/*

    Data panel data structures



    Sergei L. Kosakovsky Pond, August 2000.

*/



#ifndef _HYDATAPANEL_

#define _HYDATAPANEL_

//#pragma once

#include "HYTableWindow.h"

#include "HYSequencePanel.h"

#include "HYCanvas.h"

#include "HYPullDown.h"

#include "site.h"

#include "preferences.h"



#define  HY_DATAPANEL_NUCDATA       1

#define  HY_DATAPANEL_PROTDATA      2

#define  HY_DATAPANEL_BINARYDATA    4

#define  HY_DATAPANEL_MENU_ID       8001

#define  HY_DATAPANEL_HMENU_ID      210

#define  HY_DATAPANEL_ICON_ID       6000

#define  HY_DATAPANEL_DEF_COLORS    10

#define  HY_DATAPANEL_THERM_HSPACE  8

#define  HY_DATAPANEL_THERM_VSPACE  4

#define  HY_DATAPANEL_THERMWIDTH    19

#define  HY_DATAPANEL_CONSENSUS     1

#define  HY_DATAPANEL_RATECLASS     2

#define  HY_DATAPANEL_TRANSLATION   4

#define  HY_DATAPANEL_REFERENCE     8

#define  HY_DATAPANEL_CODEMASK      0xFFFF0000

#define  HY_DATAPANEL_OFFSETMASK    0x03

#define  HY_DATAPANEL_REVMASK       0x04

#define  HY_DATAPANEL_MODELID       0x0000FFFF

#define  HY_DATAPANEL_OPTIONS       0x000F0000

#define  HY_DATAPANEL_FREQS         0x00F00000

#define  HY_DATAPANEL_RATES         0xFF000000

#define  HY_DATAPANEL_MODEL_GLOBAL  0x01

#define  HY_DATAPANEL_MODEL_GLOBALG 0x02

#define  HY_DATAPANEL_MODEL_EFVEST  0x04

#define  HY_DATAPANEL_MODEL_MODELS  0x08

#define  HY_DATAPANEL_MODEL_NOLOCAL 0x10

#define  HY_DATAPANEL_MAX_CLASSES   32



#define  HY_WINDOW_KIND_DATAPANEL   3





//__________________________________________________________________



class _HYDataPanel: public _HYTWindow

{



public:



    _HYDataPanel(_String&,_String&);

    // data set name, tree string (or tree ident)



    virtual ~_HYDataPanel();



    virtual bool        ProcessEvent            (_HYEvent*);

    virtual bool        ProcessGEvent           (_HYEvent*);

    bool        BuildDataPanel          (void);

    void        BuildThermometer        (_HYRect* = nil);

    void        BuildMarksPane          (void);

    void        SetDataSetReference     (_String&, _SimpleList* speciesFilter = nil);

    _HYRect     ComputeNavRect          (void);

    void        SetNavRectCenter        (long,long);

    virtual void        Activate                (void);

    virtual void        Update                  (Ptr);

    virtual void        Paint                   (Ptr);

    virtual bool        ConfirmClose            (void);

    virtual bool        _ProcessMenuSelection   (long);

    virtual void        _SetMenuBar             (void);

    virtual void        _UnsetMenuBar           (void);

    virtual void        SetFilePath             (_String& s) {

        filePath=s;

    }

    virtual void        SetSavePath             (_String& s) {

        savePath=s;

    }

    virtual void        CreatePartition         (_SimpleList&,char unitLength = 1, bool jump2New = false, _String* prefix = nil);

    virtual void        KillPartition           (long);

    virtual _String     SavePartition           (long, _String* = nil, _DataSet* = nil);

    virtual void        SplitPartition          (long,long);

    virtual void        JoinPartitions          (long,long);

    virtual void        SubtractPartitions      (long,long);

    virtual void        CombPartition           (long);

    virtual void        InterleavePartitions    (long,long,bool,bool);

    virtual void        JoinSpeciesDisplay      (void);

    virtual void        OmitSelectedSpecies     (void);

    virtual void        RestoreOmittedSequence  (long);

    void        ProcessContextualPopUpMain

    (long,long);

    void        ProcessContextualPopUpAux

    (long,long);



    virtual char        WindowKind      (void)  {

        return HY_WINDOW_KIND_DATAPANEL;

    }

    bool        EditPartitionProperties (long, char changeType = 0);

    void        NavBarDblClick          (long);

    bool        SaveDataPanel           (bool saveAs = false, _String* = nil, _String* = nil, bool = true, _DataSet* = nil, bool = false);

    void        InputPartitionString    (void);

    void        RestorePartInfo         (_String*);

    void        RestorePanelSettings    (_String*);

    void        BuildLikelihoodFunction (_String* lName = nil, _SimpleList* = nil, long = -1);

    long        SpawnLikelihoodFunction (_DataSet*,_String*,_List&, _SimpleList&, _SimpleList* = nil, _SimpleList* = nil);

    long        SpawnLikelihoodFunctionNP

    (_List&, bool = false);

    void        RestoreSavedLFs         (void);

    void        ComputeLikelihoodFunction (long);

    void        SimulateDataSet         (long, bool=false);

    void        OptimizeLikelihoodFunction

    (void);



    bool        IsSelectionNonEmpty     (void);



    virtual bool        _ProcessOSEvent         (Ptr);

    void        _PrintData              (void);



    long        GetLFID                 (void) {

        return lfID;

    }

    long        GetDSID                 (void) {

        return dataSetID;

    }

    _String*        LFSnapshot              (void);

    bool        LFRestore               (long);



    void        SetHypothesis           (_String*,bool);

    long        GetHypothesis           (bool);

    long        FindLFState             (_String);

    _String     GetLFStateName          (long);

    _String*        GetLFStateString        (long);



    virtual bool        IsSaveEnabled           (void) {

        return true;

    }

    virtual bool        IsPrintEnabled          (void) {

        return true;

    }

    virtual void        SetFont                 (_HYFont&);

    void        HandleFontChange        (void);



    void        SetLockState            (bool r) {

        cantDeleteLF = r;

    }

    _DataSet*   GenerateOrderedDataSet  (void);

    void        RefreshCategoryVars     (void);

    void        FindFunction            (void);

    bool        DependOnTree              (long tid) {

        return treeVarReferences.Find (tid) >=0;

    }







private:

    void        GetMutantCount            (long,_SimpleList&,_String* = nil);

    void        ShowMutantCount           (long);

    void        CleanupSingletons         (long);

    void        InferTopologies           (bool = false);

    bool        GenerateGoodPartitions    (_SimpleList&);

    bool        CanSplit                  (void);

    bool        CanJoin                   (long,long);

    bool        CanJoinSpecies            (_SimpleList&);

    bool        CanSubtract               (long,long);

    bool        CanComb                   (long);

    bool        CanInterleave             (long, long);

    void        SelectPartition           (void);

    void        InvertSelection           (void);

    void        PartitionPropsMenu        (void);

    void        _UpdateSelectionChoices   (bool);

    void        _UpdatePartitionOperations(_SimpleList*);

    virtual void        _OmitSelectedSpecies      (_SimpleList&);

    virtual void        _RestoreOmittedSequence   (long);

    void        _CopySelectionToClipboard (void);

    void        UpdatePartitionOperations (void);

    void        UpdateSelDepPartitionOperations

    (void);

    void        _PaintThermRect           (bool update = false);

    void        _UpdateLFMenu             (void);

    void        _VerifyInferMenu          (void);

    void        BuildDataPartitions       (void);

    void        AddPartition              (_DataSetFilter*);

    long        FindUnusedColor           (void);

    void        MarkSites                 (_SimpleList&,long);

    void        UnmarkSites               (_SimpleList&,long);

    void        CorrectSites              (long);

    void        UpdateDataWrapper         (void);

    void        GenerateStatusLine        (void);

    void        GenerateTreeList          (_List&);

    void        GenerateModelList         (_List&, long);

    _String     TreeTopologyChange        (long, _String*);

    bool        DataTypeChange            (long, long);

    void        ModelChange               (_String*, long, long);

    void        ModelOptionChange         (_String*, long);

    void        ModelFreqChange           (_String*, long);

    void        ModelRateClassChange      (long,long);

    void        ShowConstantSites         (bool, bool relaxed = false, bool sequences = false);

    void        ShowDuplicateSequences    (bool);

    void        ShowCodonUsage            (long);

    void        ShowCharacterUsage        (long, bool);

    void        ShowAssociation           (long, char);

    bool        AdjustStatusLine          (long, bool force = false, long preselected = -1);

    void        ActivateInfoLines         (bool);

    void        AdjustInfoNames           (void);

    void        ExecuteProcessor          (long);

    void        UpdateConsensusSequence   (_String&, bool apply = false);

    void        UpdateTranslationString   (_String&, long,bool apply = false, long genCode = -1);

    bool        GetTranslationString      (_String&, long, char = 0,long = -1, _List * = nil);

    static void LongToPartData            (long,char&,bool&,long&);

    static long PartDataToLong            (char,bool,long);

    static void LongToModelData           (long,int&,int&,int&,int&);

    static long ModelDataToLong           (int,int,int,int);

    static void CodeTo3AA                 (_String&, long, _SimpleList*,_Parameter * = nil);

    static char CodeToAA                  (long, _SimpleList*,_Parameter * = nil);

    void        _PaintLFStatus            (void);

    bool        PurgeLF                   (bool all = true);

    bool        PurgeLFFilter             (long);

    void        DisplayParameterTable     (void);

    void        OpenGeneralBSWindow       (void);

    _String*    GetExclusionsFromCode     (long);

    _String*    GetExclusionsFromList     (_SimpleList*);

    _String*    GetExclusionsFromExcList  (_SimpleList*);

    _String*    GetMatrixFromCode         (long);

    void        IndexToCodon              (long,char*);

    void        SetCodonExclusions        (_DataSetFilter*, long, long);

    long        FindGeneticCodeByExcl     (_String*);

    void        DFExclusionsToString      (_DataSetFilter*, _String&);

    long        AddPartitionRow           (_String*,long);

    void        DeletePartitionRow        (long);

    void        ConstructDataTypeOptions  (_List&);

    void        RefreshPartRow            (_List&, long, bool = false, bool = true);

    void        GenerateModelFOptionList  (_List&, long);

    void        GenerateModelPOptionList  (_List&, long);

    void        HandleSearchAndReplace    (void);



    long        dataSetID,

                referenceSequence,

                translatedSequence,

                genCodeID,

                lfID;



    _SimpleList dataPartitions,

                partData,

                siteAssignments,

                treeVarReferences,

                modelReferences,

                partitionColors,

                overlaps,

                omittedSeqs,

                inferCacheDF;



    _List       statusLines,

                statusData,

                savedLFNames,

                inferCache,

                savedLFStates;



    _String*    dataSetName,

                filePath,

                savePath,

                topConstr;



    _HYRect     navRect,

                thermRect;



    _DataSetFilter*

    dataWrapper;



    char        dataType,

                addedLines;



    bool        tainted,

                cantDeleteLF;

    friend      class   _HYParameterTable;

    friend      class   _HYGeneralBootstrapWindow;

};





void        ReadGeneticCodes       (void);

void        NewGeneticCodeTable    (long);

void        ReadModelTemplates     (void);

void        KillLFRecordFull       (long);

long        DimensionOfGenCode     (long);



_List*      FindModelTemplate      (_String*);

long        FindModelTemplate      (_String*,long);

_List*      FindModelTemplate      (long, long);

void        SetModelMenus          (int,_HYPullDown*, _HYPullDown*);

bool        RequestDataSetReplace  (long);

bool        RequestLFDeleteOrAlter (long);

void        ReadDataPanelProcessors(void);



extern      _List   geneticCodes,

            dataPanelProcessors;





/* each item in the list is also a list of 3 items:

   table     name

   table     description

   table     definition (as a simple list of length 64)

*/



extern      _List   modelTemplates;



/* each item in the list is also a list of 3 items:

    model   name

    model   options (SimpleList: model flags, model dimension)

    model   location - path to model file

*/



#endif



//EOF
