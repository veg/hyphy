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

#include      "baseobj.h"
#include      "alignment.h"
#include      "defines.h"
#include      "batchlan.h"
#include      "likefunc.h"
#include      "bayesgraph.h"
#include      "scfg.h"
#include      "function_templates.h"
#include      "global_object_lists.h"
#include      "mersenne_twister.h"
#include      "global_things.h"
#include      "hy_string_buffer.h"

using namespace hy_global;
using namespace hyphy_global_objects;

_List       openFileHandlesBackend;
_AVLListX   openFileHandles     (&openFileHandlesBackend);

//____________________________________________________________________________________
/* various helper functions */

void _CheckExpressionForCorrectness (_Formula& parsed_expression, _String const& exp, _ExecutionList& program, long desired_type = HY_ANY_OBJECT) {
    _String error_message;
    
    long parse_result = parsed_expression.ParseFormula (exp,program.nameSpacePrefix, &error_message);
    
    if (error_message.nonempty()) {
        throw (_String ("Failed to parse ") & exp.Enquote () & " with the following error: " & error_message);
    }
    if (parse_result != HY_FORMULA_EXPRESSION) {
        throw (exp.Enquote () & " did not parse to a simple expression");
    }
    if (parsed_expression.IsEmpty ()) {
        throw (exp.Enquote () & " parsed to an empty expression");
    }
    if (!(desired_type == HY_ANY_OBJECT || parsed_expression.ObjectClass() == desired_type)) {
        // TODO SLKP 20170704: ObjectClass will compute the expression with current values which may fail
        throw (exp.Enquote () & " did not evaluate to a " & FetchObjectNameFromType (desired_type));
    }
}

//____________________________________________________________________________________

_Variable* _CheckForExistingVariableByType (_String const& name, _ExecutionList& program, long desired_type = HY_ANY_OBJECT) {
    _String variable_id = AppendContainerName(name,program.nameSpacePrefix);
    _Variable * target_variable = FetchVar(LocateVarByName (variable_id));
    
    if (!target_variable) {
        throw (variable_id.Enquote() & " is not an existing variable");
    }

    if (!(desired_type == HY_ANY_OBJECT || target_variable->ObjectClass() == desired_type)) {
        throw (name.Enquote () & " is not of type " & FetchObjectNameFromType (desired_type));
    }

    return target_variable;
}

//____________________________________________________________________________________

_PMathObj   _ProcessAnArgumentByType (_String const& expression, long desired_type, _ExecutionList& program) {
    
    _Formula  parsed_expression;
    _CheckExpressionForCorrectness (parsed_expression, expression, program, desired_type);
    
    _PMathObj expression_result = parsed_expression.Compute(0,program.nameSpacePrefix);
    if (expression_result && (expression_result->ObjectClass() & desired_type)) {
        expression_result->AddAReference();
        return expression_result;
    }
    
    throw (expression.Enquote () & " did not evaluate to a " & FetchObjectNameFromType (desired_type));
    return nil;
}

//____________________________________________________________________________________

const _String _ProcessALiteralArgument (_String const& expression, _ExecutionList& program) {
    _PMathObj the_string = _ProcessAnArgumentByType (expression, STRING, program);
    
    _String result (*((_FString*)the_string)->theString);
    the_string->RemoveAReference();
    return result;
}

  //____________________________________________________________________________________

BaseRefConst    _GetHBLObjectByType (_String const&  source_name, long& type, long * object_index = nil) {
  long            object_type = type;
  BaseRefConst    source_object = _HYRetrieveBLObjectByName (source_name, object_type,object_index,false);
  
  if (source_object == nil) {
    throw (source_name.Enquote('\'') & " is not a " & _HYHBLTypeToText(type));
  }
  type = object_type;
  return source_object;
}

  //____________________________________________________________________________________

BaseRef    _GetHBLObjectByTypeMutable (_String const&  source_name, long& type, long * object_index = nil) {
  long            object_type = type;
  BaseRef         source_object = _HYRetrieveBLObjectByNameMutable (source_name, object_type,object_index,false);
  
  if (source_object == nil) {
    throw (source_name.Enquote('\'') & " is not a " & _HYHBLTypeToText(type));
  }
  type = object_type;
  return source_object;
}


//____________________________________________________________________________________

_Variable* _ElementaryCommand::_ValidateStorageVariable (_ExecutionList& program, unsigned long argument_index) const {
    _String  storage_id (program.AddNameSpaceToID(*GetIthParameter(argument_index)));
    
    _Variable * receptacle = CheckReceptacleCommandIDException (&AppendContainerName(storage_id,program.nameSpacePrefix),get_code(), true, false, &program);

    return receptacle;
}

//____________________________________________________________________________________

bool     _DefaultExceptionHandler (_Variable * receptacle, _String const& error, _ExecutionList& current_program) {
    if (receptacle) { // if receptacle is nil, then we have already handled the error
        receptacle->SetValue(new _MathObject, false);
    }
    current_program.ReportAnExecutionError (error);
    return false;
}


//____________________________________________________________________________________

_PMathObj    _EnsurePresenceOfKey    (_AssociativeList * dict, _String const& key, long desired_type) {
    _PMathObj value = dict->GetByKey (key, desired_type);
    if (!value) {
        throw (key.Enquote() & " was not a key associated with a " & FetchObjectNameFromType (desired_type) & "-typed value");
    }
    return value;
}

//____________________________________________________________________________________

hyFloat    _NumericValueFromKey     (_AssociativeList * dict, _String const& key, hyFloat default_value) {
    _PMathObj value = dict->GetByKey (key, NUMBER);
    if (!value) {
        return default_value;
    }
    return value->Compute()->Value();
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleDifferentiate(_ExecutionList& current_program){
  _Variable * receptacle = nil;
  current_program.advance ();
  
  try {
    receptacle = _ValidateStorageVariable (current_program);
    _String     expression = *GetIthParameter(1);
    
    _Formula parsed_expression;
    _CheckExpressionForCorrectness (parsed_expression, expression, current_program);
 
    long times = 1L;
    if (parameter_count() >= 4) {
      times = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),current_program.nameSpacePrefix);
      if (times <= 0L) {
        throw (GetIthParameter(3)->Enquote() & " (the number of times to differentiate) must be a non-negative integer");
      }
    }
    
    _Formula * derivative = parsed_expression.Differentiate(*GetIthParameter(2));
    for (; times>1 && derivative; times--) {
      _Formula * temp = derivative->Differentiate (*GetIthParameter(2));
      delete derivative;
      derivative = temp;
    }
    
    receptacle->SetFormula(*derivative);
    delete derivative;
  
  } catch (const _String& error) {
      return  _DefaultExceptionHandler (receptacle, error, current_program);
  }
  return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleFindRootOrIntegrate (_ExecutionList& currentProgram, bool do_integrate){
    
    _Variable * receptacle = nil;
    currentProgram.advance ();
  
    try {
        receptacle = _ValidateStorageVariable (currentProgram);
        _String     expression = *GetIthParameter(1);
        
        _Formula parsed_expression;
        _CheckExpressionForCorrectness (parsed_expression, expression, currentProgram);
         _Variable * target_variable = _CheckForExistingVariableByType (*GetIthParameter(2),currentProgram,NUMBER);
        
        if (!parsed_expression.DependsOnVariable(target_variable->GetAVariable())) {
            throw (expression & " does not depend on the variable " & target_variable->GetName()->Enquote());
        }
         
        _Formula  * derivative = do_integrate ? nil : parsed_expression.Differentiate (*target_variable->GetName(),false);
        
        hyFloat    lb = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),currentProgram.nameSpacePrefix),
                   ub = _ProcessNumericArgumentWithExceptions (*GetIthParameter(4),currentProgram.nameSpacePrefix);
        
        if (ub<=lb) {
            throw (_String ('[') & lb & ',' & ub & "] is not a valid interval");
            return false;
        }
        
        if (!do_integrate) {
            if (derivative) {
                receptacle->SetValue (new _Constant (parsed_expression.Newton (*derivative,target_variable, 0.0, lb, ub)),false);
            } else {
                receptacle->SetValue (new _Constant (parsed_expression.Brent (target_variable, lb, ub)), false);
            }
        } else {
            receptacle->SetValue (new _Constant (parsed_expression.Integral (target_variable, lb, ub, ub-lb>1e10)), false);
        }
        
        if (derivative) {
            delete derivative;
        }
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, currentProgram);
    }
    return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleExport(_ExecutionList& current_program){
  
  _Variable * receptacle = nil;
  current_program.advance();
  
  try {
    receptacle =    _ValidateStorageVariable (current_program);

    const _String source_name   = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
    long          object_type = HY_BL_MODEL | HY_BL_LIKELIHOOD_FUNCTION | HY_BL_DATASET_FILTER | HY_BL_HBL_FUNCTION,
                  object_index;
    
    BaseRef       source_object;
    try {
      source_object = _GetHBLObjectByTypeMutable (source_name, object_type, &object_index);
    } catch (const _String& ) {
      receptacle->SetValue(new _MathObject);
    }
    
    
    switch (object_type) {
      case HY_BL_LIKELIHOOD_FUNCTION: {
        _StringBuffer * serialized_object = new _StringBuffer (8192L);
        ((_LikelihoodFunction*)source_object)->SerializeLF (*serialized_object);
        receptacle->SetValue(new _FString (serialized_object));
        break;
      }
      case HY_BL_DATASET_FILTER: {
        receptacle->SetValue(new _FString (new _String ((_String*)((_DataSetFilter*)source_object)->toStr())));
        ReleaseDataFilterLock(object_index);
        break;
      }
      case HY_BL_MODEL: {
        _StringBuffer * serialized_object = new _StringBuffer (8192L);
        SerializeModel (*serialized_object,object_index,nil,true);
        receptacle->SetValue(new _FString (serialized_object));
        break;
      }
      case HY_BL_HBL_FUNCTION: {
        receptacle->SetValue(new _FString (new _String (ExportBFFunction (object_index))));
        break;
      }
    }
    
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (receptacle, error, current_program);
  }
  return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetDataInfo (_ExecutionList& current_program) {
  
    _Variable * receptacle = nil;
    current_program.advance();
  
    try {
      
        receptacle = _ValidateStorageVariable (current_program);
        const _String source_name = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        
        long            object_type = HY_BL_DATASET|HY_BL_DATASET_FILTER;
        BaseRefConst    source_object = _GetHBLObjectByType(source_name, object_type);
 
        _DataSetFilter const * filter_source  = object_type == HY_BL_DATASET_FILTER ? (_DataSetFilter const *)source_object : nil;
        _DataSet       const * dataset_source = filter_source ? nil : (_DataSet const *)source_object;
                
        switch (parameters.lLength) {
            case 2UL: { // get site->pattern map
                if (filter_source) {
                    receptacle->SetValue (new _Matrix (filter_source->duplicateMap),false);
                } else {
                    receptacle->SetValue (new _Matrix (dataset_source->DuplicateMap()),false);
                }
            }
            break;
                
            case 3UL: { // data parameters, or sequence string
                _String const argument = _ProcessALiteralArgument (*GetIthParameter(2),current_program);
                if (argument == _String ("CHARACTERS")) {
                    _List characters;
                    if (filter_source) {
                       unsigned long character_count = filter_source->GetDimension(true),
                        fd = filter_source->GetUnitLength();
                        
                        for (unsigned long idx = 0UL; idx < character_count; idx++) {
                            characters < new _String (filter_source->ConvertCodeToLetters (filter_source->CorrectCode(idx), fd));
                        }
                    } else {
                        _String alphabet_string = dataset_source->GetTT () ? dataset_source->GetTT ()->GetAlphabetString() : kEmptyString;
                        for (unsigned long idx = 0UL; idx < alphabet_string.length(); idx++) {
                            characters < new _String (alphabet_string (idx));
                        }
                    }
                    receptacle->SetValue (new _Matrix (characters), false);
                } else if (argument == _String ("PARAMETERS")) {
                    if (filter_source) {
                        _AssociativeList * parameterInfo = new _AssociativeList;
                        
                        (*parameterInfo) < (_associative_list_key_value){"ATOM_SIZE", new _Constant (filter_source->GetUnitLength())}
                        < (_associative_list_key_value){"EXCLUSIONS", new _FString  (filter_source->GetExclusions())}
                        < (_associative_list_key_value){"SITES_STRING", new _FString  ((_String*)filter_source->theOriginalOrder.ListToPartitionString())}
                        < (_associative_list_key_value){"SEQUENCES_STRING", new _FString  ((_String*)filter_source->theNodeMap.ListToPartitionString())};
                        
                        receptacle->SetValue (parameterInfo,false);
                        
                    } else {
                        throw (argument.Enquote('\'') & " is only available for DataSetFilter objects");
                    }
                } else if (argument == _String ("CONSENSUS")) { // argument == _String("PARAMETERS")
                    if (filter_source) {
                        receptacle->SetValue (new _FString (new _String(filter_source->GenerateConsensusString())), false);
                    } else {
                        _DataSetFilter temp;
                        _SimpleList l1, l2;
                        temp.SetFilter (dataset_source, 1, l1, l2, false);
                        receptacle->SetValue (new _FString (new _String(temp.GenerateConsensusString())), false);
                    }
                } else {
                    long seqID = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix);
                    
                    if (filter_source) {
                        if (seqID>=0 && seqID < filter_source->NumberSpecies()) {
                            receptacle->SetValue (new _FString (filter_source->GetSequenceCharacters(seqID)),false);
                        } else  if (seqID >= -4 && seqID <= -1) {
                            _SimpleList indices, map, counts;
                            long uniqueSequences = filter_source->FindUniqueSequences(indices, map, counts, -seqID - 1);
                            _AssociativeList * parameterInfo = new _AssociativeList;
                            parameterInfo->MStore ("UNIQUE_SEQUENCES",             new _Constant (uniqueSequences), false);
                            parameterInfo->MStore ("UNIQUE_INDICES",            new _Matrix (indices), false);
                            parameterInfo->MStore ("SEQUENCE_MAP",          new _Matrix (map), false);
                            parameterInfo->MStore ("UNIQUE_COUNTS",      new _Matrix  (counts), false);
                            receptacle->SetValue (parameterInfo,false);
                        }
                    } else { // filter_source
                        if (seqID>=0 && seqID < dataset_source->NoOfSpecies()) {
                            receptacle->SetValue (new _FString (dataset_source->GetSequenceCharacters(seqID)),false);
                        }
                    }
                } // else numeric cases
            }
            break;
                
            case 4UL : {
                if (filter_source) {
                    long seq  = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix),
                         site = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),current_program.nameSpacePrefix);
                    
                    if (site >=0 && site<filter_source->GetPatternCount()) {
                        if ( seq>=0 && seq<filter_source->NumberSpecies()) {
                            _Matrix             * res = new _Matrix (filter_source->GetDimension (true), 1, false, true);
                            
                            bool                only_the_index = hy_env::EnvVariableTrue(hy_env::get_data_info_returns_only_the_index);
                            
                            _String             character (filter_source->RetrieveState(site, seq));
                            long                theValue = filter_source->Translate2Frequencies (character, res->theData,  true);
                            
                            if (only_the_index) {
                                receptacle->SetValue (new _Constant (theValue),false);
                                DeleteObject     (res);
                            } else {
                                receptacle->SetValue (res,false);
                            }
                        } else {
                            bool count_gaps = hy_env::EnvVariableTrue(hy_env::harvest_frequencies_gap_options);
                            long filter_dimension = filter_source->GetDimension (true);
                            
                            _Matrix * accumulator = new _Matrix (filter_dimension, 1, false, true),
                                    * storage     = new _Matrix (filter_dimension, 1, false, true);
                            
                            _String *buffer = filter_source->MakeSiteBuffer();
                            
                            for (long species_index = filter_source->NumberSpecies()-1; species_index >= 0; species_index --) {
                                filter_source->RetrieveState(site,species_index,*buffer, false);
                                filter_source->Translate2Frequencies (*buffer, storage->theData,  count_gaps >= 0.5);
                                *accumulator += *storage;
                            }
                            receptacle -> SetValue (accumulator, false);
                            BatchDelete(storage, buffer);
                            
                        }
                    } else {
                        throw (_String ("Site index ") & site & " is invalid: must be in range " & "[0, " & (long)filter_source->GetPatternCount() & "]");
                    }
                } else {
                    throw ("This set of arguments is only supported for DataSetFilter objects");
                }
            }
            break;
                        
            case 5UL: {
                if (filter_source) {
                    long seq1  = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix),
                         seq2  = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),current_program.nameSpacePrefix);
                    
                    if ( seq1>=0 && seq2 >=0 && seq1< filter_source->NumberSpecies() && seq2 <filter_source->NumberSpecies()) {
                        _String const *  res_flag = GetIthParameter(4);
                        _Matrix * res;
                        
                        if (hy_env::kPairwiseCountAmbiguitiesAverage == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingAverageFrequencyAware);
                        } else if (hy_env::kPairwiseCountAmbiguitiesResolve == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolve);
                        } else if (hy_env::kPairwiseCountAmbiguitiesSkip == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingSkip);
                        } else {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolveFrequencyAware);
                        }
                        
                        receptacle->SetValue (res,false);
                    } else {
                        throw (_String (seq1).Enquote() & "," & _String (seq2).Enquote() & " is an invalid sequence pair specification.");
                    }
                } else {
                    throw ("This set of options is not supported for DataSet arguments");
                }
            }
            break;
        // switch
        }
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }

    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetInformation (_ExecutionList& current_program) {
    _Variable * receptacle = nil;
    current_program.advance();
    try {
        
        _Matrix*   result     = nil;
        receptacle = _ValidateStorageVariable (current_program);
        const _String source_name = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        
        long            object_type = HY_BL_LIKELIHOOD_FUNCTION | HY_BL_DATASET_FILTER | HY_BL_MODEL  ,
                        object_index;
        BaseRefConst    source_object = _HYRetrieveBLObjectByName (source_name, object_type,&object_index,false);
        
        
        if (source_object) {
            switch (object_type) {
                case HY_BL_LIKELIHOOD_FUNCTION: {
                    // list of ctagory variables
                    _LikelihoodFunction const * lf = (_LikelihoodFunction const *)source_object;
                    
                    _List        catVars;
                    for (unsigned long k=0UL; k<lf->GetCategoryVars().countitems(); k++) {
                        catVars << lf->GetIthCategoryVar (k)->GetName();
                    }
                    result = new _Matrix (catVars);
                }
                break;
                case HY_BL_DATASET_FILTER : {
                    result = ((_DataSetFilter const *) source_object)->GetFilterCharacters();
                }
                break;
                case HY_BL_MODEL: {
                    _SimpleList modelParms;
                    _AVLList    modelParmsA (&modelParms);
                    
                     if (IsModelOfExplicitForm (object_index)) {
                        ((_Formula const*)source_object)->ScanFForVariables(modelParmsA,false);
                    } else {
                        ((_Variable const*)source_object)->ScanForVariables(modelParmsA,false);
                        
                    }
                    _List       modelPNames;
                    
                    for (unsigned long vi=0; vi<modelParms.lLength; vi++) {
                        modelPNames << LocateVar(modelParms.lData[vi])->GetName();
                    }
                    
                    result = new _Matrix (modelPNames);
                }
                break;
            }
        } else {
            _Variable* source_object = FetchVar(LocateVarByName (source_name));
            
            if (source_object->ObjectClass()==STRING) {
                source_object    = FetchVar (LocateVarByName (_String((_String*)source_object->Compute()->toStr())));
            }
            if (source_object) {
                if (source_object->IsCategory()) {
                    _CategoryVariable * thisCV = (_CategoryVariable*)source_object;
                    thisCV->Refresh();
                    
                    _Matrix *values  = thisCV->GetValues(),
                    *weights = thisCV->GetWeights(!thisCV->IsUncorrelated());
                    
                    long size = values->GetHDim()*values->GetVDim();
                    result = new _Matrix (2,size,false,true);
                    
                    for (unsigned long k = 0UL; k< size ; k++) {
                        result->theData[k]   = values->theData[k];
                        result->theData[size+k] = weights->theData[k];
                    }
                } else {
                    if (source_object->ObjectClass()==TREE_NODE) {
                        _CalcNode* theNode = (_CalcNode*)source_object;
                        if (theNode->GetModelIndex() != HY_NO_MODEL) {
                            result = new _Matrix;
                            theNode->RecomputeMatrix (0,1,result);
                        }
                    } else {
                        if (source_object->ObjectClass() == TOPOLOGY || source_object->ObjectClass() == TREE) {
                            
                            _List* map = ((_TreeTopology*)source_object)->MapNodesToModels ();
                            _AssociativeList* return_this = new _AssociativeList();
                            
                            for (unsigned long i = 0; i < map->lLength; i++) {
                                _List * nodeInfo = (_List*) map->GetItem(i);
                                return_this->MStore(*(_String*)nodeInfo->GetItem(0), *(_String*)nodeInfo->GetItem (1));
                            }
                            result = (_Matrix*) return_this;
                            DeleteObject (map);
                        }
                    }
                    
                    if ((!result)&& source_object->ObjectClass()==NUMBER) {
                        result = new _Matrix (1,3,false,true);
                        result->theData[0]=source_object->Compute()->Value();
                        result->theData[1]=source_object->GetLowerBound();
                        result->theData[2]=source_object->GetUpperBound();
                    }
                }
            } else {
                // TODO : SLKP 20170702, check that this still works, eh?
                _String reg_exp = GetStringFromFormula (&source_name,current_program.nameSpacePrefix);
                if (reg_exp != *source_name) {
                    int errNo = 0;
                    regex_t* regex = PrepRegExp (reg_exp, errNo, true);
                    if (regex) {
                        _List       matches;
                        
                        for (AVLListXIteratorKeyValue variable_record : AVLListXIterator (&variableNames)) {
                            _String* vName = (_String*)variableNames.Retrieve (variable_record.get_index());
                            if (vName->RegExpMatch (regex).lLength) {
                                matches << vName;
                            }
                        }
                        if (matches.lLength) {
                            result = new _Matrix (matches);
                        }
                        _String::FlushRegExp (regex);
                    } else {
                        HandleApplicationError (_String::GetRegExpError (errNo));
                    }
                }
            }
        }
        if (!result) {
            result = new _Matrix (0,0,false,false);
        }
        receptacle->SetValue(result, false);
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleConstructCategoryMatrix (_ExecutionList& current_program) {
    
    static      _Trie        kRunOptions (
                                            _String ("COMPLETE"), CategoryConstructionOptions::kMatrixConditionals,
                                            _String ("WEIGHTS"), CategoryConstructionOptions::kMatrixWeights,
                                            _String ("SITE_LOG_LIKELIHOODS"), CategoryConstructionOptions::kSiteProbabilities,
                                            _String ("CLASSES"), CategoryConstructionOptions::kMatrixClasses
                                          );

    
    _Variable * receptacle = nil;
    current_program.advance();
  
    try {
        
        _Matrix*   result     = nil;
        receptacle = _ValidateStorageVariable (current_program);
        const _String source_name = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        long            object_type = HY_BL_LIKELIHOOD_FUNCTION  | HY_BL_TREE,
                        object_index;
        BaseRefConst    source_object = _GetHBLObjectByType (source_name, object_type, &object_index);
      
        switch (object_type) {
            case HY_BL_LIKELIHOOD_FUNCTION: {
                _Matrix * partition_list  = nil;
                if (parameters.countitems () > 3) { // have a restricting partition
                    partition_list = (_Matrix*)_ProcessAnArgumentByType(GetIthParameter(3), MATRIX, current_program);
                }
                _SimpleList   included_partitions;
                _LikelihoodFunction * like_func = (_LikelihoodFunction * )source_object;
                
                like_func->ProcessPartitionList(included_partitions, partition_list);
                DeleteObject (partition_list);
                
                CategoryConstructionOptions::CategoryConstructionOptions run_mode = CategoryConstructionOptions::kMatrixConditionals;
 
                if (parameters.countitems () > 2) {
                    long run_mode_long = kRunOptions.GetValueFromString(*GetIthParameter(2));
                    if (run_mode_long != kNotFound) {
                        run_mode = (CategoryConstructionOptions::CategoryConstructionOptions)run_mode_long;
                    }
                }
            
                receptacle->SetValue(like_func->ConstructCategoryMatrix(included_partitions, run_mode ,true, receptacle->GetName()), false);
            }
            break;
            case HY_BL_TREE: {
                _TheTree  *source_tree       = (_TheTree*)source_object;
                
                long    which_partition   = 0L,
                        linked_likelihood_id = source_tree->IsLinkedToALF (which_partition);
                
                if (linked_likelihood_id >= 0) {
                    _LikelihoodFunction * linked_lf            = (_LikelihoodFunction*) likeFuncList (linked_likelihood_id);
                    const _DataSetFilter      * filter             = linked_lf->GetIthFilter (which_partition);
                    linked_lf->PrepareToCompute();
                    linked_lf->Compute         ();
                    long patterns                         = filter->GetPatternCount();
                    
                    _Matrix             *conditional_matrix     = new _Matrix   (2*patterns*(source_tree->GetLeafCount()
                                                                                 + source_tree->GetINodeCount()) * source_tree->categoryCount,
                                                                     source_tree->GetCodeBase(),
                                                                     false, true);
                    
                    _List               leaf_names,
                                        internal_names;
                    
                    _TreeIterator ti (source_tree, _HY_TREE_TRAVERSAL_POSTORDER);
                    
                    while (_CalcNode * iterator = ti.Next()) {
                        if (ti.IsAtLeaf()) {
                            leaf_names.AppendNewInstance (new _String(iterator->ContextFreeName()));
                        } else {
                            internal_names.AppendNewInstance (new _String(iterator->ContextFreeName()));
                        }
                    }
                    
                    leaf_names << internal_names;
                    
                    for (unsigned long site = 0UL; site < patterns; site ++) {
                        source_tree->RecoverNodeSupportStates (filter,site,*conditional_matrix);
                    }
                    
                    linked_lf->DoneComputing   ();
                    receptacle->SetValue (
                        &(*(new _AssociativeList)
                          < _associative_list_key_value ({"Nodes", new _Matrix (leaf_names)})
                          < _associative_list_key_value ({"Values", conditional_matrix})),
                        false
                    );
                }
            }
        }
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    
    return true;
    
}
//____________________________________________________________________________________

bool      _ElementaryCommand::HandleAlignSequences(_ExecutionList& current_program) {
    _Variable * receptacle = nil;
    
    static    const   _String   kCharacterMap           ("SEQ_ALIGN_CHARACTER_MAP"),
                                kScoreMatrix            ("SEQ_ALIGN_SCORE_MATRIX"),
                                kGapChar                ("SEQ_ALIGN_GAP_CHARACTER"),
                                kGapOpen                ("SEQ_ALIGN_GAP_OPEN"),
                                kGapExtend              ("SEQ_ALIGN_GAP_EXTEND"),
                                kGapOpen2               ("SEQ_ALIGN_GAP_OPEN2"),
                                kGapExtend2             ("SEQ_ALIGN_GAP_EXTEND2"),
                                kFrameShift             ("SEQ_ALIGN_FRAMESHIFT"),
                                kGapLocal               ("SEQ_ALIGN_NO_TP"),
                                kAffineGaps             ("SEQ_ALIGN_AFFINE"),
                                kCodonAlign             ("SEQ_ALIGN_CODON_ALIGN"),
                                kLinearSpace            ("SEQ_ALIGN_LINEAR_SPACE"),
                                kScoreMatrixCodon3x1    ("SEQ_ALIGN_PARTIAL_3x1_SCORES"),
                                kScoreMatrixCodon3x2    ("SEQ_ALIGN_PARTIAL_3x2_SCORES"),
                                kScoreMatrixCodon3x4    ("SEQ_ALIGN_PARTIAL_3x4_SCORES"),
                                kScoreMatrixCodon3x5    ("SEQ_ALIGN_PARTIAL_3x5_SCORES"),
                                kLocalAlignment         ("SEQ_ALIGN_LOCAL_ALIGNMENT");


    current_program.advance();
  
    try {
        _Matrix   * result     = nil;
        receptacle = _ValidateStorageVariable (current_program);
        _Matrix   * input_seqs = (_Matrix   *)_ProcessAnArgumentByType(*GetIthParameter(1), MATRIX, current_program);
        
        unsigned long        input_seq_count = input_seqs->GetSize ();
        
        auto    string_validator = [] (long row, long col, _Formula* cell) -> bool {
            if (cell) {
                if (cell->ObjectClass() != STRING) {
                    throw (_String(" Matrix entry (") & row & "," & col & ") did not evaluate to a string");
                }
                return true;
            }
            throw (_String("Empty matrix entry (") & row & "," & col & ")");
            return false;
        };

        if (! (input_seqs->IsAStringMatrix() && (input_seqs->is_row() || input_seqs->is_column()) && input_seq_count >= 2 && input_seqs->ValidateFormulaEntries (string_validator))) {
            throw (GetIthParameter(1)->Enquote() & " did not evaluate to a dense string vector with ≥2 entries");
        }
        
        
        _AssociativeList* alignment_options  = (_AssociativeList   *)_ProcessAnArgumentByType(*GetIthParameter(2), ASSOCIATIVE_LIST, current_program);
        _FString        * char_vector        = (_FString*)          _EnsurePresenceOfKey (alignment_options, kCharacterMap, STRING);
        
        unsigned          long     char_count = 0UL;
        long              character_map_to_integers [256] = {-1L};
        
        for (unsigned long cc = 0UL; cc < char_vector->theString->length(); cc++) {
            unsigned char this_char = char_vector->theString->get_uchar(cc);
            if (character_map_to_integers [this_char]>=0) {
                throw (_String ("Duplicate character ") & _String ((char)this_char).Enquote('\'') & " in " & kCharacterMap);
            } else {
                character_map_to_integers [this_char] = cc;
                char_count ++;
            }
        }
        if (char_count == 0) {
            throw ("Null alphabet supplied");
        }
        
        bool        do_local       = _NumericValueFromKey (alignment_options,kGapLocal,0.0) > 0.5,
                    do_affine      = _NumericValueFromKey (alignment_options,kAffineGaps,0.0) > 0.5,
                    do_linear      = _NumericValueFromKey (alignment_options,kLinearSpace,1.0) > 0.5,
                    do_codon       = _NumericValueFromKey (alignment_options,kCodonAlign,0.0) > 0.5,
                    do_full_local  = do_codon && _NumericValueFromKey (alignment_options,kLocalAlignment,0.0) > 0.5;
        
        
        long        codon_count        = char_count * char_count * char_count,
                    expected_dimension = do_codon ? codon_count + 1UL : char_count + 1UL;
    
        _Matrix *   score_matrix = (_Matrix*)_EnsurePresenceOfKey(alignment_options, kScoreMatrix, MATRIX);
        
        if (!score_matrix->check_dimension (expected_dimension, expected_dimension)) {
            throw (_String ("The dimension of the scoring matrix ") & kScoreMatrix.Enquote('(',')') & " was not the expected dimension: " & expected_dimension & 'x' & expected_dimension);
        }
        
        score_matrix = (_Matrix*)score_matrix->ComputeNumeric();
        score_matrix->CheckIfSparseEnough(true); // force to not be sparse
        
        _Matrix     *codon3x5 = nil,
                    *codon3x4 = nil,
                    *codon3x2 = nil,
                    *codon3x1 = nil;
        
        if (do_codon) {
            
             unsigned long expected_columns [4] = {codon_count * 10UL,codon_count * 4UL,char_count * char_count * 3UL,char_count * 3UL};
            _String const * keys [4]            = {&kScoreMatrixCodon3x5, &kScoreMatrixCodon3x4, &kScoreMatrixCodon3x2, &kScoreMatrixCodon3x1};
            _Matrix **      targets [4]         = {&codon3x5, &codon3x4, &codon3x2, &codon3x1};
            
            for (int i = 0; i < 4; i++) {
                (*targets[i]) = (_Matrix*)_EnsurePresenceOfKey(alignment_options, *keys[i], MATRIX);
                if (!(*targets[i])->check_dimension (expected_dimension, expected_columns[i])) {
                    throw (_String ("The dimension of the scoring matrix ") & keys[i]->Enquote('(',')') & " was not the expected dimension: " & expected_dimension & 'x' & (long)expected_columns[i]);
                }
                (*targets[i]) = (_Matrix*)(*targets[i])->ComputeNumeric();
                (*targets[i]) -> CheckIfSparseEnough(true);
            }
            
            for (unsigned long i = 0UL; i < 256UL; i++) {
            // this maps all undefined characters to '?' essentially
                if (character_map_to_integers[i] < 0) {
                    character_map_to_integers[i] = -codon_count - 1;
                }
            }
        }
        
        char        gap_character = '-';
        if (_FString    *gap_c = (_FString*)alignment_options->GetByKey (kGapChar, STRING)) {
            if (gap_c->theString->length () != 1UL) {
                throw (_String ("Invalid gap character specification ") &  *gap_c->theString);
            }
            gap_character = gap_c->theString->char_at(0UL);
        }
        
        hyFloat  gap_open       = _NumericValueFromKey (alignment_options, kGapOpen,15.),
                 gap_open2      = _NumericValueFromKey (alignment_options, kGapOpen2,gap_open),
                 gap_extend     = _NumericValueFromKey (alignment_options, kGapExtend,1.),
                 gap_extend2    = _NumericValueFromKey (alignment_options, kGapExtend2,gap_extend),
                 gap_frameshift = _NumericValueFromKey (alignment_options, kFrameShift,50.);
        
        _StringBuffer settings_report (256L);
        
        settings_report << "\n\tGap character               : " << gap_character
                        << "\n\tGap open cost [reference]   : " << gap_open
                        << "\n\tGap open cost [query]       : " << gap_open2
                        << "\n\tGap extend cost [reference] : " << gap_extend
                        << "\n\tGap extend cost [query]     : " << gap_extend2
                        << "\n\tCodon frameshift cost       : " << gap_frameshift
                        << "\n\tIgnore terminal gaps        : " << (do_local?"Yes":"No")
                        << "\n\tPerform local alignment     : " << (do_full_local?"Yes":"No");
    
        if (do_codon) {
            settings_report << "\n\tUse codon alignment with frameshift routines";
            do_linear = false;
        }
        
        _AssociativeList *aligned_strings = new _AssociativeList;
        _String * reference_sequence = ((_FString*)input_seqs->GetFormula(0,0)->Compute())->theString;
        
        
        for (unsigned long index2 = 1UL; index2 < input_seq_count; index2++) {
            _String * sequence2 = ((_FString*)input_seqs->GetFormula(0,index2)->Compute())->theString;
            _AssociativeList * pairwise_alignment = new _AssociativeList;
            hyFloat    score = 0.0;
            if (do_linear) {
                unsigned long   size_allocation = sequence2->length()+1UL;
                
                _Matrix         *buffers[6];
                
                ArrayForEach(buffers, 6, [=] (_Matrix* m, unsigned long) -> _Matrix* {
                    return new _Matrix (size_allocation,1,false,true);
                });

                char          *alignment_route = new char[2*(size_allocation)] {0};
                
                _SimpleList ops (reference_sequence->length()+2UL,-2,0);
                ops.lData[reference_sequence->length()+1] = sequence2->length();
                ops.lData[0]               = -1;
                
                score = LinearSpaceAlign(reference_sequence,sequence2,character_map_to_integers,score_matrix,
                                         gap_open,gap_extend,gap_open2,gap_extend2,
                                         do_local,do_affine,ops,score,0,
                                         reference_sequence->length(),0,sequence2->length(),buffers,0,alignment_route);
                
                delete[]    alignment_route;
                
                _StringBuffer     *result1 = new _StringBuffer (reference_sequence->length() + 1UL),
                                  *result2 = new _StringBuffer (size_allocation);
                
                long             last_column     = ops.lData[ops.lLength-1];
                
                for (long position = (long)reference_sequence->length() - 1L; position>=0; position--) {
                    long current_column     = ops.lData[position+1];
                    
                    if (current_column<0) {
                        if (current_column == -2 /*|| (current_column == -3 && last_column == string2->sLength)*/) {
                            current_column = last_column;
                        } else if (current_column == -3) {
                            // find the next matched char or a -1
                            long    p   = position, s2p;
                            while ((ops.lData[p+1]) < -1) {
                                p--;
                            }
                            
                            s2p = ops.lData[p+1];
                            
                            for (long j = last_column-1; j>s2p;) {
                                (*result1) << gap_character;
                                (*result2) << sequence2->char_at(j--);
                            }
                            
                            last_column     = s2p+1;
                            
                            for (; position>p; position--) {
                                (*result2) << gap_character;
                                (*result1) << reference_sequence->char_at(position);
                            }
                            position ++;
                            continue;
                        } else {
                            for (last_column--; last_column >=0L; last_column--) {
                                (*result1) << gap_character;
                                (*result2) << sequence2->char_at (last_column);
                            }
                            while (position>=0) {
                                (*result1) << reference_sequence->char_at (position--);
                                (*result2) << gap_character;
                            }
                            break;
                        }
                    }
                    
                    if (current_column == last_column) { // insert in sequence 2
                        (*result1) << reference_sequence->char_at (position);
                        (*result2) << gap_character;
                    } else {
                        last_column--;
                        
                        for (; last_column > current_column; last_column--) { // insert in column 1
                            (*result2) << sequence2->char_at (last_column);
                            (*result1) << gap_character;
                        }
                        (*result1) << reference_sequence->char_at (position);
                        (*result2) << sequence2->char_at (current_column);
                    }
                    //printf ("%s\n%s\n", result1->sData, result2->sData);
                }
                
                for (last_column--; last_column >=0; last_column--) {
                    (*result1) << gap_character;
                    (*result2) << sequence2->char_at(last_column);
                }
                
                result1->Flip ();
                result2->Flip ();
                pairwise_alignment->MStore ("1", new _FString(result1), false);
                pairwise_alignment->MStore ("2", new _FString(result2), false);
            } else { // not linear
                char * str1r = nil,
                     * str2r = nil;
                
                score = AlignStrings (reference_sequence->get_str(),
                                      sequence2->get_str(),
                                      str1r,
                                      str2r,
                                      character_map_to_integers,
                                      score_matrix->fastIndex(),
                                      score_matrix->GetVDim(),
                                      gap_character,
                                      gap_open,
                                      gap_extend,
                                      gap_open2,
                                      gap_extend2,
                                      gap_frameshift,
                                      do_local,
                                      do_affine,
                                      do_codon,
                                      char_count,
                                      do_codon ? codon3x5->fastIndex() : nil,
                                      do_codon ? codon3x4->fastIndex() : nil,
                                      do_codon ? codon3x2->fastIndex() : nil,
                                      do_codon ? codon3x1->fastIndex() : nil,
                                      do_full_local);
                
                if ( str1r && str2r ) {
                    pairwise_alignment->MStore ("1", new _FString (new _String( str1r )), false);
                    pairwise_alignment->MStore ("2", new _FString (new _String( str2r )), false);
                    delete [] str1r;
                    delete [] str2r;
                } else {
                    throw ( "Internal Error in AlignStrings" );
                }
                pairwise_alignment->MStore ("0", new _Constant (score), false);
                aligned_strings->MStore (_String((long)index2-1L), pairwise_alignment, false);
            }
        }
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    return true;
}
//____________________________________________________________________________________

bool      _ElementaryCommand::HandleHarvestFrequencies (_ExecutionList& current_program) {
    
    _Variable * receptacle = nil;
    current_program.advance();
  
    try {
        _Matrix   * result     = nil;
        
        receptacle = _ValidateStorageVariable (current_program);

        long       object_type = HY_BL_DATASET|HY_BL_DATASET_FILTER;
        BaseRefConst    source_object = _GetHBLObjectByType(*GetIthParameter(1), object_type);
      
        long      unit      = _ProcessNumericArgumentWithExceptions(*GetIthParameter(2),current_program.nameSpacePrefix),
                  atom      = _ProcessNumericArgumentWithExceptions(*GetIthParameter(3),current_program.nameSpacePrefix);
        
        bool      position_specific = _ProcessNumericArgumentWithExceptions(*GetIthParameter(3),current_program.nameSpacePrefix) > 0.5,
                  include_gaps       = hy_env::EnvVariableTrue(hy_env::harvest_frequencies_gap_options);
        
        switch (object_type) {
            case HY_BL_DATASET: {
                _String vertical_partition      (parameters.countitems () > 5 ? *GetIthParameter(5) : kEmptyString),
                        horizontal_partition    (parameters.countitems () > 6 ? *GetIthParameter(6) : kEmptyString);
                
                _DataSet const * dataset = (_DataSet const*)source_object;
                _SimpleList     processed_sequence_partition, processed_site_partition;
                dataset->ProcessPartition (horizontal_partition,processed_sequence_partition,false);
                dataset->ProcessPartition (vertical_partition,processed_site_partition,true);

                receptacle->SetValue (dataset->HarvestFrequencies(unit,atom,position_specific,processed_sequence_partition, processed_site_partition,include_gaps), false);
            }
            break;
            case HY_BL_DATASET_FILTER: {
                receptacle->SetValue (((_DataSetFilter const*)source_object)->HarvestFrequencies(unit,atom,position_specific,include_gaps), false);
            }
            break;
        }
        
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleOptimizeCovarianceMatrix (_ExecutionList& current_program, bool do_optimize) {
    _Variable * receptacle = nil;
    current_program.advance();
  
    try {
        receptacle = _ValidateStorageVariable (current_program);
        
        long       object_type = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM|HY_BL_HBL_FUNCTION;
        _String    optimize_me = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        
        _LikelihoodFunction*   source_object = (_LikelihoodFunction*)_HYRetrieveBLObjectByNameMutable (AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix), object_type,nil,do_optimize==false);
        
        if (!source_object) { // Custom function (expression based)
            source_object = new _CustomFunction (optimize_me, current_program.nameSpacePrefix);
        }
        
        if (do_optimize) {
            receptacle -> SetValue(source_object->Optimize(),false);
        } else {
            _PMathObj     covariance_parameters = hy_env::EnvVariableGet(hy_env::covariance_parameter, ASSOCIATIVE_LIST|STRING);
            _SimpleList   *restrictor = nil;
            switch (object_type) {
                case HY_BL_LIKELIHOOD_FUNCTION:
                case HY_BL_SCFG: {
                    if (covariance_parameters) { // only consider some variables
                        _SimpleList variable_ids;
                        if (covariance_parameters->ObjectClass () == ASSOCIATIVE_LIST) {
                            // a list of variables stored as keys in an associative array
                            restrictor = new _SimpleList;
                            _List*  restricted_variables = ((_AssociativeList*)covariance_parameters)->GetKeys();
                            for (unsigned long iid = 0; iid < restricted_variables->lLength; iid++) {
                                variable_ids << LocateVarByName (current_program.AddNameSpaceToID(*(_String*)(*restricted_variables)(iid)));
                            }
                            DeleteObject (restricted_variables);
                        } else { // STRING
                            variable_ids << LocateVarByName (current_program.AddNameSpaceToID(*((_FString*)covariance_parameters)->theString));
                        }
                        if (!variable_ids.empty()) {
                            restrictor = new _SimpleList ();

                            for (unsigned long var_index = 0UL; var_index < variable_ids.lLength; var_index++) {
                                // TODO SLKP 20170706: this is a very inefficient linear search over and over again (in a large array)
                                long vID = source_object->GetIndependentVars().Find(variable_ids (var_index));
                                if (vID >= 0) (*restrictor) << vID;
                            }
                            
                            if (restrictor->empty() == 0) {
                                DeleteAndZeroObject (restrictor);
                            }
                        }
                    }
                }
                break;
                case HY_BL_BGM: {
                    _Matrix * bgm_cov = (_Matrix*)source_object->CovarianceMatrix(nil);
                    if (bgm_cov) {
                        receptacle->SetValue(bgm_cov,false);
                    } // TODO SLKP 20170706: handle the case when null is returned (why would that happen?); warn the user.
                }
                break;
            }
            
            receptacle->SetValue(source_object->CovarianceMatrix(restrictor),false);
            DeleteObject (restrictor);
        }
        
        if (object_type == HY_BL_NOT_DEFINED) {
            DeleteObject (source_object);    // delete the custom function object
        }
        
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleReplicateConstraint (_ExecutionList& current_program) {
    // TODO SLKP 20170706 this needs to be reimplemented; legacy code is ugly and buggy
    _Variable * receptacle = nil;
    
    try {
        
        current_program.currentCommand++;
        _String  const  constraint_pattern = _ProcessALiteralArgument (*GetIthParameter(0), current_program);
 /*
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
  HandleApplicationError (*(_String*)parameters(0)&" has no 'this' references in call to ReplicateConstraint!");
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
  HandleApplicationError (_String("Invalid reference to ") & thisS & " in the constraint specification");
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
  HandleApplicationError (_String("Unused ") & ind1 & "-th reference variable: " & *(_String*)parameters(ind1));
  return ;
  }
  
  _String namespaced = chain.AddNameSpaceToID(*(_String*)parameters(ind1));
  ind2 = LocateVarByName (namespaced);
  if (ind2<0) {
  HandleApplicationError(namespaced & " is undefined in call to ReplicateConstraint.");
  return  ;
  }
  
  _Variable* thisNode = FetchVar (ind2);
  if (thisNode->ObjectClass()==TREE_NODE) {
  thisArgs<< (long)((_CalcNode*)thisNode)->LocateMeInTree();
  } else if (thisNode->ObjectClass()==TREE) {
  thisArgs<< (long)&((_TheTree*)thisNode)->GetRoot();
  } else {
  HandleApplicationError (*(_String*)parameters(ind1) & " is neither a tree nor a tree node in call to ReplicateConstraint.");
  return ;
  }
  }
  
  // now with this list ready we can recurse down the tree and produce the contsraints
  if (RecurseDownTheTree(thisArgs, parameters, theConstraints, parts, thisIndex)) {
  if (theConstraints.lLength) {
  ReportWarning  (_String("\nReplicateConstraint generated the following contsraints:"));
  hyFloat      doDeferSet;
  checkParameter (deferConstrainAssignment,doDeferSet,0.0);
  bool            applyNow = CheckEqual(doDeferSet,0.0);
  _String         *constraintAccumulator = new _String(128L,true);
  
  if (applyNow) {
  deferSetFormula = new _SimpleList;
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
*/
  
        // first find the pattern
  
  
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
  
    return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleComputeLFFunction (_ExecutionList& current_program) {
  
  const static _String kLFStartCompute ("LF_START_COMPUTE"),
                       kLFDoneCompute  ("LF_DONE_COMPUTE");
  
  current_program.advance();
  _Variable * receptacle = nil;
  
  try {
    
    
    _String    const op_kind = * GetIthParameter(1UL);
    
    long       object_type = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM;
    _LikelihoodFunction*    source_object = (_LikelihoodFunction*)_GetHBLObjectByType(AppendContainerName (*GetIthParameter(0UL), current_program.nameSpacePrefix),object_type);
    
    if (op_kind == kLFStartCompute) {
      source_object->PrepareToCompute(true);
    } else if (op_kind == kLFDoneCompute) {
      source_object->DoneComputing (true);
    } else {
      if (!source_object->HasBeenSetup()) {
        throw (_String("Please call LFCompute (, ") & *GetIthParameter (0UL)& kLFStartCompute & ") before evaluating the likelihood function");
      } else {
        receptacle = _ValidateStorageVariable (current_program, 2);
        receptacle->SetValue (new _Constant (source_object->Compute()), false);
      }
    }
    
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (receptacle, error, current_program);
  }
  
  return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleUseModel (_ExecutionList& current_program) {
  const static _String kUseNoModel ("USE_NO_MODEL");
  current_program.advance();

  try {
    
    _String    raw_model_name = *GetIthParameter(0UL),
               source_name  = AppendContainerName (raw_model_name, current_program.nameSpacePrefix);

    long       object_type_request = HY_BL_MODEL,
               object_type = object_type_request,
               model_index = HY_NO_MODEL;
    
    _Matrix*    source_model = (_Matrix*)_HYRetrieveBLObjectByNameMutable (source_name, object_type,&model_index,false);
    
    if (!source_model && raw_model_name != kUseNoModel) {
        throw (source_name.Enquote() & " does not refer to a valid defined substitution model and is not " & useNoModel);
    }
    
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  
  return true;
}

//____________________________________________________________________________________
bool      _ElementaryCommand::HandleRequireVersion(_ExecutionList& current_program){
  current_program.advance();
  
  try {
    _String requested_version = _ProcessALiteralArgument (*GetIthParameter(0UL),current_program);
    if (kHyPhyVersion.to_float() < requested_version.to_float()) {
      throw (_String ("Current script requires at least version ")& requested_version &" of HyPhy. Please download an updated version from http://www.hyphy.org or github.com/veg/hyphy and try again.");
    }
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
}


  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleDeleteObject(_ExecutionList& current_program){
  current_program.advance();
  
  for (unsigned long i = 0UL; i < parameter_count(); i++) {
    long       requested_type = HY_BL_LIKELIHOOD_FUNCTION,
               object_index   = kNotFound;
    BaseRef    source_object = _HYRetrieveBLObjectByNameMutable (AppendContainerName(*GetIthParameter(i),current_program.nameSpacePrefix), requested_type,&object_index,false);
    
    if  (source_object) {
      KillLFRecord (object_index,true);
    } else {
      ReportWarning(GetIthParameter(i)->Enquote() & " is not a supported agrument type for " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(get_code()));
    }
  }
  return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleClearConstraints(_ExecutionList& current_program){
  current_program.advance();

  for (unsigned long i = 0UL; i< parameter_count(); i++) {
    
    _String source_name  = AppendContainerName (*(_String*)parameters(i), current_program.nameSpacePrefix);
    _Variable *clear_me = (_Variable*) FetchVar (LocateVarByName (source_name));
    if (clear_me) { // variable exists
      clear_me->ClearConstraints();
    } else {
      ReportWarning(GetIthParameter(i)->Enquote() & " is not an existing variable in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(get_code()));
    }
  }
  return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetURL(_ExecutionList& current_program){
  const static _String   save_to_file_action  ("SAVE_TO_FILE");
  
  current_program.advance();
  _Variable * receptacle = nil;
  
  try {
    _String url = _ProcessALiteralArgument (*GetIthParameter(1UL),current_program),
            *action = GetIthParameter(2UL, false);
    
    if (!action) { // store URL contents in a variable
      receptacle = _ValidateStorageVariable (current_program);
      if (Get_a_URL(url)) {
        receptacle->SetValue(new _FString (url,false),false);
      } else {
        throw (_String ("Could not fetch ") & url.Enquote());
      }
    } else {
      if (*action == save_to_file_action) {
        _String file_name = _ProcessALiteralArgument (*GetIthParameter(1UL),current_program);
        if (!ProcessFileName (file_name, true,true,(hyPointer)current_program.nameSpacePrefix),false,&current_program) {
          return false;
        }
        if (!Get_a_URL(url, &file_name)) {
          throw (_String ("Could not fetch ") & url.Enquote());
        }
      } else {
        throw (_String("Unknown action ") & action->Enquote());
      }
    }
    
    
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleAssert (_ExecutionList& current_program) {
  current_program.advance();

  try {
    _Formula parsed_expression;
    _CheckExpressionForCorrectness (parsed_expression, *GetIthParameter(0UL), current_program, NUMBER);
    if (CheckEqual (parsed_expression.Compute()->Value (), 0.0)) { // assertion failed
      bool soft_assertions = hy_env::EnvVariableTrue(hy_env::assertion_behavior);
      _String assertion_feedback;
      _String * custom_error_message = GetIthParameter(1UL,false);
      if (custom_error_message) {
        assertion_feedback = _ProcessALiteralArgument(*custom_error_message, current_program);
      } else {
        assertion_feedback = _String("Assertion ") & GetIthParameter(0UL)->Enquote() & " failed.";
      }
      if (soft_assertions) {
        StringToConsole (assertion_feedback);
        NLToConsole();
        current_program.GoToLastInstruction ();
      } else {
        current_program.ReportAnExecutionError (assertion_feedback);
      }
    }
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
}


  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleSelectTemplateModel (_ExecutionList& current_program) {
  static _String last_model_used;
  
  current_program.advance();
  try {
    _String source_name = *GetIthParameter(0UL);
    if (source_name == hy_env::use_last_model) {
      if (last_model_used.nonempty()) {
        PushFilePath (last_model_used);
      } else {
        throw ( hy_env::use_last_model &" cannot be used before any models have been defined.");
      }
    } else {
      ReadModelList();
      
      long            object_type = HY_BL_DATASET|HY_BL_DATASET_FILTER;
      _DataSetFilter const *    source_filter = (_DataSetFilter const*)_GetHBLObjectByType(source_name, object_type);

      
      _String             data_type;
      unsigned long       unit_length      = source_filter->GetUnitLength();
      
      _TranslationTable const*  filter_table = source_filter->GetData()->GetTT();
      
      if (unit_length==1UL) {
        if (filter_table->IsStandardNucleotide()) {
          data_type = "nucleotide";
        } else if (filter_table->IsStandardAA()) {
          data_type = "aminoacid";
        }
      } else {
        if (filter_table->IsStandardNucleotide()) {
          if (unit_length==3UL) {
            data_type = "codon";
          } else {
            if (unit_length==2UL)
              data_type = "dinucleotide";
          }
        }
      }
      
      if (data_type.empty()) {
        throw (source_name.Enquote () & " contains non-standard data and template models can't be selected on it");
      }

      _SimpleList matching_models;

      for (unsigned long model_index = 0; model_index < templateModelList.lLength; model_index++) {
        _List *model_components = (_List*)templateModelList(model_index);
        
        if (data_type == *(_String*)model_components->GetItem(3)) {
          _String * dim = (_String*)model_components->GetItem(2);
          if (*dim== _String("*")|| source_filter->GetDimension() == dim->to_long()) {
            matching_models << model_index;
          }
        }
      }
      
      if (matching_models.empty()) {
        throw (source_name.Enquote () & " could not be matched with any template models");
      }

      long model_id = kNotFound;
      
      if (current_program.stdinRedirect) {
        _String const option = current_program.FetchFromStdinRedirect ();
        
        model_id = matching_models.FindOnCondition( [&] (long index) -> bool {
          return option == * (_String*) templateModelList.GetItem (index,0);
        });
        
        
        if (model_id == kNotFound) {
          throw (option.Enquote() & " is not a valid model (with input redirect)");
          return false;
        }
      } else {
        #ifdef __HEADLESS__
          throw ("Unhandled standard input interaction in SelectTemplateModel for headless HyPhy");
        #endif
        
        
        
        for (int i = 0; i < kMaxDialogPrompts; i++) {
          printf ("\n\n               +--------------------------+\n");
          printf (    "               | Select a standard model. |\n");
          printf (    "               +--------------------------+\n\n\n");
          
          for (model_id = 0; model_id<matching_models.lLength; model_id++) {
            printf ("\n\t(%s):%s",((_String*)templateModelList.GetItem((matching_models(model_id),0)))->get_str(),
                                  ((_String*)templateModelList.GetItem((matching_models(model_id),1)))->get_str());
          }
          printf ("\n\n Please type in the abbreviation for the model you want to use:");
          _String const user_choice = StringFromConsole();
 
          model_id = matching_models.FindOnCondition( [&] (long index) -> bool {
            return user_choice.EqualIgnoringCase(*(_String*) templateModelList.GetItem (index,0));
          });
          
          if (model_id != kNotFound) {
            break;
          }
        }
        
        if (model_id == kNotFound) {
          throw ("Dialog did not return a valid choice after maximum allowed number of tries");
          return false;
        }

      }
      
      _String  model_file = GetStandardDirectory (HY_HBL_DIRECTORY_TEMPLATE_MODELS) & *(_String*)templateModelList.GetItem((matching_models(model_id),4));

      _ExecutionList       std_model;
      PushFilePath        (model_file, false);
      ReadBatchFile       (model_file,std_model);
      PopFilePath         ();
      lastModelUsed       = model_file;
      std_model.Execute (&current_program);
      
    }
  } catch (const _String& error) {
      return  _DefaultExceptionHandler (nil, error, current_program);
  }
  
  return true;
  
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleMolecularClock(_ExecutionList& current_program){
  current_program.advance();
  try {
    _CalcNode * apply_clock_here = (_CalcNode *)_CheckForExistingVariableByType (*GetIthParameter(0), current_program, TREE | TREE_NODE);
    _TheTree  * parent_tree;
    
    _String   clock_base;
    
    if (apply_clock_here->ObjectClass() == TREE_NODE) {
      parent_tree = (_TheTree*)((_VariableContainer*)apply_clock_here)->GetTheParent();
      if (!parent_tree) {
        throw (_String("Internal error - orphaned tree node ") & apply_clock_here->GetName ()->Enquote());
      }
      clock_base = apply_clock_here->GetName() -> Cut (parent_tree->GetName()->length() + 1, kStringEnd);
    } else {
      parent_tree = (_TheTree*)apply_clock_here;
    }

    parent_tree->MolecularClock(clock_base,parameters);

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
}


//____________________________________________________________________________________

bool      _ElementaryCommand::HandleSetParameter (_ExecutionList& current_program) {

  current_program.advance();

  try {
    _String object_to_change = *GetIthParameter(0UL);
    
    /* handle special cases */
    if (object_to_change == hy_env::random_seed) {
      hy_env::EnvVariableSet(hy_env::random_seed, new _Constant (hy_random_seed = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix) ), false);
      return true;
    }

    if (object_to_change == hy_env::defer_constrain_assignment) {
      bool defer_status = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix);
      if (defer_status) {
        deferSetFormula = new _SimpleList;
      } else if (deferSetFormula) {
        FinishDeferredSF ();
      }
      return true;
    }

    if (object_to_change == hy_env::random_seed) {
      current_program.errorHandlingMode = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix);
      return true;
 
    }

    if (object_to_change == hy_env::status_bar_update_string) {
      SetStatusLineUser (_ProcessALiteralArgument (*GetIthParameter(1), current_program));
      return true;
      
    }

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
 
  
  _String *currentArgument = (_String*)parameters(0),
  nmspc           = AppendContainerName(*currentArgument,currentProgram.nameSpacePrefix),
  errMsg,
  result;
  
  

  long objectIndex,
  typeFlag    = HY_BL_ANY;
  
  BaseRef theObject      = _HYRetrieveBLObjectByNameMutable (nmspc, typeFlag, &objectIndex);
  
  switch (typeFlag)
  {
    case HY_BL_BGM: { // BGM Branch
      currentArgument = (_String*)parameters(1);
      
      _BayesianGraphicalModel * lkf = (_BayesianGraphicalModel *) theObject;
        // set data matrix
      if (currentArgument->Equal (&bgmData)) {
        _Matrix     * dataMx = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);
        if (dataMx) {
          long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();
          
          if (dataMx->GetVDim() == num_nodes) {
            ((_BayesianGraphicalModel *)lkf)->SetDataMatrix ((_Matrix *) dataMx);
          } else {
            currentProgram.ReportAnExecutionError (_String("Data matrix columns (") & dataMx->GetVDim() & " ) does not match number of nodes in graph (" & num_nodes & ")");
            return false;
          }
        } else {
          return false;
        }
        
      }
      
        // restore node score cache
      else if (currentArgument->Equal (&bgmScores)) {
        _AssociativeList * cacheAVL = (_AssociativeList *)FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
        if (cacheAVL) {
          ((_BayesianGraphicalModel *)lkf)->ImportCache (cacheAVL);
        } else {
          return false;
        }
      }
      
        // set structure to user-specified adjacency matrix
      else if (currentArgument->Equal (&bgmGraph)) {
        _Matrix     * graphMx   = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);
        
        if (graphMx) {
          long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();
          
          if (graphMx->GetHDim() == num_nodes && graphMx->GetVDim() == num_nodes) {
            ((_BayesianGraphicalModel *)lkf)->SetStructure ((_Matrix *) graphMx->makeDynamic());
          } else {
            currentProgram.ReportAnExecutionError("Dimension of graph does not match current graph");
            return false;
          }
        } else {
          return false;
        }
      }
      
        // set constraint matrix
      else if (currentArgument->Equal (&bgmConstraintMx)) {
        _Matrix     * constraintMx  = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);
        
        if (constraintMx) {
          long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();
          
          if (constraintMx->GetHDim() == num_nodes && constraintMx->GetVDim() == num_nodes) {
            ((_BayesianGraphicalModel *)lkf)->SetConstraints ((_Matrix *) constraintMx->makeDynamic());
          } else {
            currentProgram.ReportAnExecutionError ("Dimensions of constraint matrix do not match current graph");
            return false;
          }
        } else {
          return false;
        }
      }
      
        // set node order
      else if (currentArgument->Equal (&bgmNodeOrder)) {
        _Matrix     * orderMx   = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);
        
        if (orderMx) {
            // UNDER DEVELOPMENT  April 17, 2008 afyp
          long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();
          
          _SimpleList     * orderList = new _SimpleList();
          
          orderList->Populate (num_nodes, 0, 0);
          
          if (orderMx->GetVDim() == num_nodes) {
            for (long i = 0; i < num_nodes; i++) {
              orderList->lData[i] = (long) ((*orderMx) (0, i));
            }
            
            ((_BayesianGraphicalModel *)lkf)->SetNodeOrder ( (_SimpleList *) orderList->makeDynamic() );
          } else {
            currentProgram.ReportAnExecutionError ("Length of order vector doesn't match number of nodes in graph");
            return false;
          }
        } else {
          return false;
        }
      }
      
      
        // set network parameters
      else if (currentArgument->Equal (&bgmParameters)) {
        _AssociativeList* inAVL = (_AssociativeList*)FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
        if (inAVL) {
          ((_BayesianGraphicalModel *)lkf)->ImportCache (inAVL);
        } else {
          return false;
        }
      }
      
      
        // anything else
      else {
        currentProgram.ReportAnExecutionError (*currentArgument & " is not a valid BGM parameter");
        return false;
      }
    } // end BGM
      break;
      
    case HY_BL_SCFG:
    case HY_BL_LIKELIHOOD_FUNCTION:
    {
      currentArgument = (_String*)parameters(1);
      if (typeFlag == HY_BL_SCFG && currentArgument->Equal (&scfgCorpus)) {
        ((Scfg*)theObject)->SetStringCorpus ((_String*)parameters(2));
      } else {
        _LikelihoodFunction * lkf = (_LikelihoodFunction *) theObject;
        currentArgument = (_String*)parameters(1);
        long g = ProcessNumericArgument(currentArgument,currentProgram.nameSpacePrefix);
        if (lkf->GetIndependentVars().Map (g) < 0L) {
          currentProgram.ReportAnExecutionError (*currentArgument & " (=" & g & ") is not a valid parameter index");
          return false;
        }
        currentArgument = (_String*)parameters(2);
        lkf->SetIthIndependent (g,ProcessNumericArgument(currentArgument,currentProgram.nameSpacePrefix));
      }
    }
      break;
        // end SCFG and LF
      
    case HY_BL_DATASET:
    case HY_BL_DATASET_FILTER: {
      _DataSet * ds = nil;
      
      
      long sequence_index  = ProcessNumericArgument ((_String*)parameters(1),currentProgram.nameSpacePrefix);
      if (typeFlag == HY_BL_DATASET) {
        ds = (_DataSet*) theObject;
      }
      else {
        _DataSetFilter *dsf = (_DataSetFilter*)theObject;
        ds = dsf->GetData ();
        sequence_index = dsf->theNodeMap.Map (sequence_index);
      }
      
      if (typeFlag == HY_BL_DATASET_FILTER) {
        ReleaseDataFilterLock(objectIndex);
      }
      
      
      _String * sequence_name = new _String(ProcessLiteralArgument ((_String*)parameters(2),currentProgram.nameSpacePrefix));
      
      if (! ds->SetSequenceName (sequence_index, sequence_name)) {
        delete sequence_name;
        currentProgram.ReportAnExecutionError (*((_String*)parameters(1)) & " (=" & sequence_index & ") is not a valid sequence index");
        return false;
        
      }
    } // end data set and data set filter
      break;
        // Dataset and Datasetfilter
      
    default:
        // check to see if this is a calcnode
      _CalcNode* treeNode = (_CalcNode*)FetchObjectFromVariableByType(&nmspc, TREE_NODE);
      if (treeNode) {
        if (*((_String*)parameters(1)) == _String("MODEL")) {
          _String modelName = AppendContainerName(*((_String*)parameters(2)),currentProgram.nameSpacePrefix);
          long modelType = HY_BL_MODEL, modelIndex;
          BaseRef modelObject      = _HYRetrieveBLObjectByNameMutable (modelName, modelType, &modelIndex, true);
          if (modelObject) {
            _VariableContainer * parentTree = treeNode->ParentTree();
            if (!parentTree) {
              currentProgram.ReportAnExecutionError (*((_String*)parameters(0)) & " is an orphaned tree node (the parent tree has been deleted)");
              return false;
            }
            long pID, lfID = ((_TheTree*)parentTree->Compute())->IsLinkedToALF(pID);
            if (lfID>=0){
              currentProgram.ReportAnExecutionError ((*parentTree->GetName()) & " is linked to a likelihood function (" & *GetObjectNameByType (HY_BL_LIKELIHOOD_FUNCTION, lfID) &") and cannot be modified ");
              return false;
            }
            
            treeNode->ReplaceModel (modelName, parentTree);
            break;
          }
          else {
            currentProgram.ReportAnExecutionError (*((_String*)parameters(2)) & " does not appear to be a valid model name");
            return false;
          }
        } else {
          currentProgram.ReportAnExecutionError (*((_String*)parameters(1)) & " is not a supported parameter type for a tree node argument");
          return false;
        }
      }
      
      currentProgram.ReportAnExecutionError (*currentArgument & " is not a valid likelihood function/data set filter/tree topology/tree node");
      return false;
      
      
  } // end cases
  return true;
}


//____________________________________________________________________________________
//____________________________________________________________________________________












//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetString (_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
  
    _String  errMsg,
             *result = nil;

    long    f,
            sID,
            sID2 = -1;

    _Variable * theReceptacle = CheckReceptacleCommandID (&AppendContainerName(*((_String*)parameters(0)),currentProgram.nameSpacePrefix),HY_HBL_COMMAND_GET_STRING, true, false, &currentProgram);

    if (!theReceptacle) {
        return false;
    }

    sID = ProcessNumericArgument ((_String*)parameters(2),currentProgram.nameSpacePrefix);
    if (parameters.lLength>3) {
        sID2 = ProcessNumericArgument ((_String*)parameters(3),currentProgram.nameSpacePrefix);
    }

    f = _HY_GetStringGlobalTypes.Find((_String*)parameters(1));
    if (f >=0 ) {
        f = _HY_GetStringGlobalTypes.GetXtra (f);
    }
  
    switch (f) {

        case HY_BL_LIKELIHOOD_FUNCTION: // LikelihoodFunction
        case HY_BL_DATASET:
        case HY_BL_DATASET_FILTER:
        case HY_BL_SCFG:
        case HY_BL_BGM: {
            result = (_String*)GetObjectNameByType(f,sID);
            if (result) {
                result = (_String*) result->makeDynamic();
				//ReportWarning(_String((const char*)"In HandleGetString(): ") & result);
            }
            break;
        }
           
        case HY_BL_HBL_FUNCTION: // UserFunction
            result = (_String*)GetObjectNameByType(HY_BL_HBL_FUNCTION,sID);
            if (result) {
                _AssociativeList * resAVL = new _AssociativeList;
                resAVL->MStore ("ID", new _FString (*result), false);
                resAVL->MStore ("Arguments", new _Matrix(GetBFFunctionArgumentList(sID)), false);
                theReceptacle->SetValue (resAVL,false);
                return true;
            } 
            break;
            
        case HY_BL_TREE: { // Tree
            // 20110608 SLKP: REFACTOR into a separate function
            // I am sure this is used somewhere else (perhaps for other types)
            result = FetchMathObjectNameOfTypeByIndex (TREE, sID);
            if (result) {
                result = (_String*) result->makeDynamic();
            }
            break;
        }

        default: { // everything else...
            // decide what kind of object current argument represents
          
          
          
            _String *currentArgument = (_String*)parameters(1),
                    nmspaced       = AppendContainerName(*currentArgument,currentProgram.nameSpacePrefix);
            long    typeFlag       = HY_BL_ANY,
                    index          = -1;
                    
            BaseRefConst theObject      = _HYRetrieveBLObjectByName (nmspaced, typeFlag, &index);

            if (theObject) {
                switch (typeFlag) {
                case HY_BL_DATASET: {
                    _DataSet const* dataSetObject = (_DataSet const*)theObject;
                    if (sID>=0 && sID<dataSetObject->NoOfSpecies()) {
                        result = (_String*)(dataSetObject->GetNames().GetItem(sID))->makeDynamic();
                    } else {
                        if (sID < 0) {
                            theReceptacle->SetValue (new _Matrix (dataSetObject->GetNames()), false);
                            return true;
                        }
                    }
                    break;
                }
                case HY_BL_DATASET_FILTER: {
                    _DataSetFilter const* dataSetFilterObject = (_DataSetFilter const*)theObject;

                    if (sID >=0 && sID<dataSetFilterObject->NumberSpecies()) {
                        result = (_String*)dataSetFilterObject->GetData()->GetNames().GetItem(dataSetFilterObject->theNodeMap.Element(sID))->makeDynamic();
                    } else {
                        if (sID < 0) {
                            _List filterSeqNames;
                            _List const * originalNames = &dataSetFilterObject->GetData()->GetNames();
                          
                            for (long seqID=0; seqID<dataSetFilterObject->NumberSpecies(); seqID++) {
                                filterSeqNames << originalNames->GetItem (dataSetFilterObject->theNodeMap.Element (seqID));
                            }
                            theReceptacle->SetValue (new _Matrix (filterSeqNames), false);
                            return true;
                        }
                    }
                    break;
                }
                case HY_BL_BGM: {
                    //ReportWarning(_String("In HandleGetString() for case HY_BL_BGM"));
                    _BayesianGraphicalModel * this_bgm      = (_BayesianGraphicalModel *) theObject;

                    switch (sID) {
                        case HY_HBL_GET_STRING_BGM_SCORE: {   // return associative list containing node score cache
                            _AssociativeList        * export_alist  = new _AssociativeList;

                            if (this_bgm -> ExportCache (export_alist)) {
                                theReceptacle -> SetValue (export_alist, false);
                                return true;
                            } else {
                                DeleteObject (export_alist);
                                errMsg = _String ("Failed to export node score cache for BGM '") & nmspaced & "'";
                            }

                            break;
                        }
                        case HY_HBL_GET_STRING_BGM_SERIALIZE: {   // return associative list with network structure and parameters
                            result = new _String (1024L, true);
                            this_bgm -> SerializeBGM (*result);
                            result->Finalize();
                            theReceptacle->SetValue (new _FString (result),false);
                            return true;
                        }
                        default: {
                            errMsg = _String ("Unrecognized index ") & sID & " for a BGM object";
                        }
                    }
                }
                case HY_BL_LIKELIHOOD_FUNCTION:
                case HY_BL_SCFG: {
                
                    _LikelihoodFunction *lf = (_LikelihoodFunction*)theObject;
                    if (sID>=0) {
                        if (sID<lf->GetIndependentVars().lLength) {
                            result = (_String*)(LocateVar(lf->GetIndependentVars().lData[sID])->GetName())->makeDynamic();
                        } else {
                            if (sID<lf->GetIndependentVars().lLength+lf->GetDependentVars().lLength) {
                                result = (_String*)(LocateVar(lf->GetDependentVars().lData[sID-lf->GetIndependentVars().lLength])->GetName())->makeDynamic();
                            }
                        }
                    } else {
                        _AssociativeList* resList = lf->CollectLFAttributes ();
                        if (typeFlag == HY_BL_SCFG) {
                            ((Scfg*)lf)->AddSCFGInfo (resList);
                        }
                        theReceptacle->SetValue (resList,false);
                        return true;
                    }
                    break;
                }
                
                case HY_BL_MODEL: {
                    if (sID>=0) {
                        // check to make see if the 
                       if (sID2 < 0) { // get the sID's parameter name
                            _SimpleList     modelP;
                            _AVLList        modelPA (&modelP);
                            ScanModelForVariables (index, modelPA, false, -1, false);
                            modelPA.ReorderList();
                            if (sID<modelP.lLength) {
                                result = (_String*)LocateVar(modelP.lData[sID])->GetName()->makeDynamic();
                            } 

                        } else { // get the formula for cell (sID, sID2)
                            if (!IsModelOfExplicitForm (index)) {
                                _Variable*      theMx = (_Variable*)theObject;
                                _Formula * cellFla = ((_Matrix*)theMx->GetValue())->GetFormula (sID,sID2);
                                if (cellFla) {
                                    result = new _String((_String*)cellFla->toStr());
                                }
                            }
                        }
                        
                    } else {
                        _Variable   * tV, * tV2;
                        bool         mByF;
                        RetrieveModelComponents (index, tV, tV2, mByF);

                        if (tV) {
                            if (sID == -1) { // branch length expression
                                result = ((_Matrix*)tV->GetValue())->BranchLengthExpression((_Matrix*)tV2->GetValue(),mByF);
                            } else
                                /*
                                    returns an AVL with keys
                                    "RATE_MATRIX" - the ID of the rate matrix
                                    "EQ_FREQS"    - the ID of eq. freq. vector
                                    "MULT_BY_FREQ" - a 0/1 flag to determine which format the matrix is in.
                                */
                            {
                                _AssociativeList * resList = new _AssociativeList;
                                resList->MStore ("RATE_MATRIX",new _FString(*tV->GetName()),false);
                                resList->MStore ("EQ_FREQS",new _FString(*tV2->GetName()),false);
                                resList->MStore ("MULT_BY_FREQ",new _Constant (mByF),false);
                                theReceptacle->SetValue (resList,false);
                                return true;
                            }
                        }
                    }
                    break;
                }
                case HY_BL_HBL_FUNCTION: {
                    _AssociativeList * resAVL = new _AssociativeList;
                    resAVL->MStore ("ID", new _FString (*GetObjectNameByType (HY_BL_HBL_FUNCTION, index, false)), false);
                    resAVL->MStore ("Arguments", new _Matrix(GetBFFunctionArgumentList(index)), false);
                    resAVL->MStore ("Body", new _FString (GetBFFunctionBody(index).sourceText,false),false);
                    theReceptacle->SetValue (resAVL,false);
                    return true;
                }
            } // end of "switch" 
        }
        else {
            if (currentArgument->Equal(&versionString)) {
                if (sID > 1.5)
#ifdef __HEADLESS__
                    result = new _String(_String ("Library version ") & kHyPhyVersion);
#else
#ifdef __MAC__
                    result = new _String(_String("Macintosh ") & kHyPhyVersion);
#else
#ifdef __WINDOZE__
                    result = new _String(_String("Windows ") & kHyPhyVersion);
#else
                    result = new _String(_String("Source ") & kHyPhyVersion);
#endif
#endif
#endif
                    else if (sID > 0.5) {
                        result = new _String(GetVersionString());
                    } else {
                        result = new _String(kHyPhyVersion);
                    }
                } else if (currentArgument->Equal(&timeStamp)) {
                    result = new _String(GetTimeStamp(sID < 0.5));
                } else if (currentArgument->Equal(&listLoadedLibraries)) {
                  theReceptacle->SetValue (new _Matrix (loadedLibraryPaths.Keys()));
                  return true;
                } else {
                  _Variable* theVar = FetchVar(LocateVarByName (nmspaced));
                  if (theVar) {
                    if (theVar->IsIndependent() && sID != -3) {
                      result = (_String*)theVar->toStr();
                    } else {
                      if (sID == -1){
                        
                        _SimpleList vL;
                        _AVLList    vAVL (&vL);
                        theVar->ScanForVariables (vAVL, true);
                        vAVL.ReorderList();
                        _AssociativeList   * resL = new _AssociativeList;
                        _List splitVars;
                        SplitVariableIDsIntoLocalAndGlobal (vL, splitVars);
                        InsertVarIDsInList (resL, "Global", *(_SimpleList*)splitVars(0));
                        InsertVarIDsInList (resL, "Local",   *(_SimpleList*)splitVars(1));
                        
                        theReceptacle->SetValue (resL,false);
                        return true;
                      }
                      
                      else {  // formula string
                        
                        if (sID == -3) {
                          _String local, global;
                          _SimpleList var_index;
                          var_index << theVar->GetAVariable ();
                          if (theVar->IsIndependent()) {
                            //printf ("ExportIndVariables\n");
                            ExportIndVariables (global, local, &var_index);
                          } else {
                            //printf ("ExportDepVariables\n");
                            ExportDepVariables(global, local, &var_index);
                          }
                          result = new _String (128L, true);
                          (*result) << global << local << '\n';
                          result->Finalize();
                        } else {
                          _Matrix * formula_matrix = (sID2 >= 0 && theVar->ObjectClass() == MATRIX) ? (_Matrix*)theVar->GetValue () : nil;
                          if (formula_matrix) {
                           _Formula* cell = formula_matrix->GetFormula(sID, sID2);
                           if (cell) {
                              result = (_String*) cell->toStr();
                            }
                          } else {
                            result = (_String*)theVar->GetFormulaString ();
                          }
                        }
                      }
                    }
                  } else {
                    errMsg = _String ("'") & *currentArgument & "' is not an allowed argument type ";
                  }
                }
            }
        }
    }

    if (errMsg.sLength) {
        currentProgram.ReportAnExecutionError (errMsg); 
        DeleteObject (result);
        result = nil;    
    }
    
    if (result) {
        theReceptacle->SetValue (new _FString (result),false);
        return true;
    }

    theReceptacle->SetValue (new _MathObject(), false);
    

    return false;
}






//____________________________________________________________________________________

bool      _ElementaryCommand::HandleFprintf (_ExecutionList& currentProgram)
{
    currentProgram.currentCommand++;
    _String* targetName = (_String*)parameters(0),
             fnm;
    
    bool     doClose                 = true,
             print_to_stdout         = false;
    
    FILE*   dest = nil;
    
    try {
        bool    skipFilePathEval        = false;
        
        if (targetName->Equal(&stdoutDestination)) {
            _FString * redirect = (_FString*)FetchObjectFromVariableByType (&blFprintfRedirect, STRING);
            if (redirect && redirect->theString->sLength) {
                if (redirect->theString->Equal (&blFprintfDevNull)) {
                    return true; // "print" to /dev/null
                } else {
                    skipFilePathEval = true;
                    targetName       = redirect->theString;
                } 
            }
            else {
                print_to_stdout = true;
            }
        }
        
        checkParameter (printDigitsSpec,printDigits,0L);
        
        if (!print_to_stdout) {
            fnm = *targetName;
            if (fnm.Equal(&messageLogDestination)) {
                if ((dest = hy_message_log_file) == nil) {
                    return true; // requested print to MESSAGE_LOG, but it does not exist
                                 // (e.g. remote MPI nodes, or running from a read only location
                }
            } else {
                if (skipFilePathEval == false && !fnm.IsALiteralArgument()) {
                    fnm = GetStringFromFormula (&fnm,currentProgram.nameSpacePrefix);
                }
                
                if (!fnm.ProcessFileName(true,false,(hyPointer)currentProgram.nameSpacePrefix, false, &currentProgram)) {
                    return false;
                }
                
                
                long k  = openFileHandles.Find (&fnm);
                
                doClose = k<0;
                
                if (!doClose) {
                    dest = (FILE*)openFileHandles.GetXtra (k);
                } else {
                    if ((dest = doFileOpen (fnm.getStr(), "a")) == nil)
                        throw  (_String  ("Could not create/open output file at path '") & fnm & "'.");
                }
            }
        }
        
        for (long i = 1; i<parameters.lLength; i++) {
            _String    *varname = (_String*)parameters(i);
            
            BaseRef    thePrintObject   =   nil;
            _Formula   f;
            
            if (varname->Equal(&clearFile)) {
                if (!print_to_stdout && dest) {
                    fclose (dest);
                    dest = doFileOpen (fnm.getStr(), "w");
                    long k = openFileHandles.Find (&fnm);
                    if (k>=0) {
                        openFileHandles.SetXtra(k, (long)dest);
                    }
                }
            } else if (varname->Equal(&keepFileOpen) && !print_to_stdout) {
                if (openFileHandles.Find (&fnm) < 0) {
                    openFileHandles.Insert (fnm.makeDynamic(), (long)dest);
                }
                doClose = false;
            } else if (varname->Equal(&closeFile) && !print_to_stdout) {
                openFileHandles.Delete (&fnm, true);
                doClose = true;
            } else if (varname->Equal(&systemVariableDump)) {
                thePrintObject=&variableNames;
            } else if (varname->Equal(&selfDump)) {
                thePrintObject=&currentProgram;
            } else {
                // check for possible string reference
                
                _String    temp    = ProcessStringArgument (varname),
                           nmspace;
                
                if (temp.sLength > 0) {
                    nmspace = AppendContainerName(temp,currentProgram.nameSpacePrefix);
                    if (nmspace.IsValidIdentifier(fIDAllowCompound)) {
                        thePrintObject = FetchObjectFromVariableByType (&nmspace,HY_ANY_OBJECT);
                    }
                } else {
                    nmspace = AppendContainerName(*varname,currentProgram.nameSpacePrefix);
                }
                
                
                if (thePrintObject == nil) {
                    long typeFlag = HY_BL_ANY,
                         index;
                  
                    thePrintObject = _HYRetrieveBLObjectByNameMutable (nmspace, typeFlag, &index);
                    
                    if (!thePrintObject) {
                        _String argCopy = *varname,
                                errMsg;

                        _FormulaParsingContext fpc (&errMsg, currentProgram.nameSpacePrefix);
                        if (Parse (&f,argCopy, fpc, nil) == HY_FORMULA_EXPRESSION) {
                            thePrintObject = f.Compute(0,currentProgram.nameSpacePrefix);
                        } else {
                            if (errMsg.sLength)
                                throw (errMsg);
                            else
                                throw (_String ("Argument ") & i & " is not a simple expression");
                        }
                    } else {
                      if (typeFlag == HY_BL_DATASET_FILTER) {
                        ReleaseDataFilterLock(index);
                      }
                    }
                }
            }
            
            if (thePrintObject) {
                if (!print_to_stdout) {
                    thePrintObject->toFileStr (dest);
                } else {
                    _String outS ((_String*)thePrintObject->toStr());
                    StringToConsole (outS);
                }
            }
        }
    }
    catch (_String errMsg) {
        currentProgram.ReportAnExecutionError (errMsg);
    }
    
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    if (print_to_stdout) {
        yieldCPUTime();
    }
#endif
    if (dest && dest!=hy_message_log_file && doClose) {
        fclose (dest);
    }
    
    return !currentProgram.IsErrorState();
}


