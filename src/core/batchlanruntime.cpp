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
#include      "global_object_lists.h"
#include      "mersenne_twister.h"
#include      "global_things.h"
#include      "hy_string_buffer.h"
#include      "associative_list.h"
#include      "tree_iterator.h"

#include      "function_templates.h"


#ifndef __HYPHY_NO_SQLITE__
  #include "sqlite3.h"
#endif

using namespace hy_global;
using namespace hyphy_global_objects;

#include <ctype.h>
//____________________________________________________________________________________
/* various helper functions */

const _String   _ElementaryCommand::ExtractStatementAssignment (_String const& source, long& end_at, const bool validate, const bool exceptions, const long offset) {
    
    _String id;
    
    try {
        long id_start = source.FirstNonSpaceFollowingSpace(offset);
        
        end_at = id_start != kNotFound ? source.FindTerminator(id_start, '=') : kNotFound;
        
        if (id_start == kNotFound || end_at == kNotFound) {
            throw _String ("Missing an ID in 'Type <ID> = statement'");
        }
        
        id = source.Cut (id_start, end_at - 1);
        if (validate) {
            if (!id.IsValidIdentifier(fIDAllowCompound)) {
                throw id.Enquote() & "is not a valid storage variable identifier";
            }
        }
        
        end_at ++;
        
    } catch (const _String& err) {
        if (exceptions) {
            throw err;
        }
        id.Clear();
        end_at = kNotFound;
    }
    
    return id;
    
}

//____________________________________________________________________________________

const _String   _ElementaryCommand::ProcessProcedureCall (_String const& source,_String& procedure, _List& pieces) {
    long op_start;
    _String id = ExtractStatementAssignment (source, op_start, false);
    
    long paren_start = op_start,
    paren_end  = source.ExtractEnclosedExpression(paren_start, '(', ')', fExtractRespectQuote | fExtractRespectEscape);
    
    if (paren_end == kNotFound) {
        throw _String ("Missing () enclosed argument list");
    }
    
    procedure = source.Cut (op_start,paren_start-1L);
    pieces < new _String (id);
    ExtractConditions (source,paren_start+1,pieces,',');
    
    return id;
}

//____________________________________________________________________________________

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
    if (!(desired_type == HY_ANY_OBJECT || parsed_expression.ObjectClass() & desired_type)) {
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

    if (!(desired_type == HY_ANY_OBJECT || target_variable->ObjectClass() & desired_type)) {
        throw (name.Enquote () & " is not of type " & FetchObjectNameFromType (desired_type));
    }

    return target_variable;
}

//____________________________________________________________________________________

HBLObjectRef   _ProcessAnArgumentByType (_String const& expression, long desired_type, _ExecutionList& program, _List* reference_manager) {
    // The return value needs to managed by the caller

    /* first see if this is a simple expression of the form 'variable_id' */

    HBLObjectRef simple_var = FetchObjectFromVariableByType (&AppendContainerName (expression, program.nameSpacePrefix), desired_type);
    if (simple_var) {
      if (reference_manager)
        *reference_manager << simple_var;
      else
        simple_var->AddAReference();
      return simple_var;
    }

    _Formula  parsed_expression;
    _CheckExpressionForCorrectness (parsed_expression, expression, program, desired_type);

    HBLObjectRef expression_result = parsed_expression.Compute(0,program.nameSpacePrefix);
    if (expression_result && (expression_result->ObjectClass() & desired_type)) {
        if (reference_manager)
          *reference_manager << expression_result;
        else
          expression_result->AddAReference();
        return expression_result;
    }

    throw (expression.Enquote () & " did not evaluate to a " & FetchObjectNameFromType (desired_type));
    return nil;
}

//____________________________________________________________________________________

const _String _ProcessALiteralArgument (_String const& expression, _ExecutionList& program) {
    HBLObjectRef the_string = _ProcessAnArgumentByType (expression, STRING, program, nil);

    _String result (((_FString*)the_string)->get_str());
    DeleteObject (the_string);
    return result;
}

  //____________________________________________________________________________________

BaseRefConst    _GetHBLObjectByType (_String const&  source_name, long& type, long * object_index, _ExecutionList* current_program) {
  long            object_type = type;
  BaseRefConst    source_object = _HYRetrieveBLObjectByName (source_name, object_type,object_index,false, true, current_program);

  if (source_object == nil) {
    throw (source_name.Enquote('\'') & " is not a " & _HYHBLTypeToText(type));
  }
  type = object_type;
  return source_object;
}

  //____________________________________________________________________________________

BaseRef    _GetHBLObjectByTypeMutable (_String const&  source_name, long& type, long * object_index = nil, bool do_literal_lookup = true) {
  long            object_type = type;
  BaseRef         source_object = _HYRetrieveBLObjectByNameMutable (source_name, object_type,object_index,false, do_literal_lookup);

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
        receptacle->SetValue(new _MathObject, false, true, NULL);
    }
    current_program.ReportAnExecutionError (error);
    return false;
}


//____________________________________________________________________________________

HBLObjectRef    _EnsurePresenceOfKey    (_AssociativeList * dict, _String const& key, long desired_type) {
    HBLObjectRef value = dict->GetByKey (key, desired_type);
    if (!value) {
        throw (key.Enquote() & " was not a key associated with a " & FetchObjectNameFromType (desired_type) & "-typed value");
    }
    return value;
}

//____________________________________________________________________________________

hyFloat    _NumericValueFromKey     (_AssociativeList * dict, _String const& key, hyFloat default_value) {
    HBLObjectRef value = dict->GetByKey (key, NUMBER);
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
    if (parameter_count() >= 4UL) {
      times = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),current_program.nameSpacePrefix);
      if (times <= 0L) {
        throw (GetIthParameter(3UL)->Enquote() & " (the number of times to differentiate) must be a non-negative integer");
      }
    }
      
    _Variable * dx  = _ValidateStorageVariable (current_program, 2);
    _Formula * derivative = parsed_expression.Differentiate(*dx->GetName());
      
    for (; times>1 && derivative; times--) {
      _Formula * temp = derivative->Differentiate (*GetIthParameter(2));
      delete derivative;
      derivative = temp;
    }
    if (derivative) {
      receptacle->SetFormula(*derivative);
      delete derivative;
    } else {
      throw (_String ("Differentiation of ") & _String((_String*)GetIthParameter(1)->toStr()).Enquote() & " failed.");
    }

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
        _Variable * target_variable = _ValidateStorageVariable (currentProgram, 2); // create variable if it doesn't exist
        target_variable = _CheckForExistingVariableByType (*GetIthParameter(2),currentProgram,NUMBER);

        if (!parsed_expression.DependsOnVariable(target_variable->get_index()) && !do_integrate) {
            throw (expression.Enquote() & " does not depend on the variable " & target_variable->GetName()->Enquote());
        }

 
        hyFloat    lb = _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),currentProgram.nameSpacePrefix),
                   ub = _ProcessNumericArgumentWithExceptions (*GetIthParameter(4),currentProgram.nameSpacePrefix);

        if (ub<=lb) {
            throw (_String ('[') & lb & ',' & ub & "] is not a valid interval");
            return false;
        }

        _Formula  * derivative = do_integrate ? nil : parsed_expression.Differentiate (*target_variable->GetName(),false);

        if (!do_integrate) {
            if (derivative) {
                receptacle->SetValue (new _Constant (parsed_expression.Newton (*derivative,target_variable, 0.0, lb, ub)),false,true, NULL);
            } else {
                receptacle->SetValue (new _Constant (parsed_expression.Brent (target_variable, lb, ub)), false,true, NULL);
            }
        } else {
            receptacle->SetValue (new _Constant (parsed_expression.Integral (target_variable, lb, ub, ub-lb>100)), false, true, NULL);
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
      receptacle->SetValue(new _MathObject, false, true, NULL);
    }


    switch (object_type) {
      case HY_BL_LIKELIHOOD_FUNCTION: {
        _StringBuffer * serialized_object = new _StringBuffer (8192L);
        ((_LikelihoodFunction*)source_object)->SerializeLF (*serialized_object);
        receptacle->SetValue(new _FString (serialized_object), false, true, NULL);
        break;
      }
      case HY_BL_DATASET_FILTER: {
        receptacle->SetValue(new _FString (new _String ((_String*)((_DataSetFilter*)source_object)->toStr())), false, true, NULL);
        ReleaseDataFilterLock(object_index);
        break;
      }
      case HY_BL_MODEL: {
        _StringBuffer * serialized_object = new _StringBuffer (8192L);
        SerializeModel (*serialized_object,object_index,nil,true);
        receptacle->SetValue(new _FString (serialized_object), false, true, NULL);
        break;
      }
      case HY_BL_HBL_FUNCTION: {
        receptacle->SetValue(new _FString (new _String (ExportBFFunction (object_index))), false, true, NULL);
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
     static const _String kPairwiseCountAmbiguitiesResolve                ("RESOLVE_AMBIGUITIES"),
                          kPairwiseCountAmbiguitiesAverage                ("AVERAGE_AMBIGUITIES"),
                          kPairwiseCountAmbiguitiesSkip                   ("SKIP_AMBIGUITIES"),
                          kCharacters                                     ("CHARACTERS"),
                          kConsensus                                      ("CONSENSUS"),
                          kParameters                                     ("PARAMETERS"),
                          kPattern                                        ("PATTERN"),
                          kSite                                           ("SITE");


    _Variable * receptacle = nil;
    current_program.advance();

    try {

        receptacle = _ValidateStorageVariable (current_program);
        const _String source_name = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);

        long            object_type = HY_BL_DATASET|HY_BL_DATASET_FILTER;
        BaseRefConst    source_object = _GetHBLObjectByType(source_name, object_type, nil, &current_program);

        _DataSetFilter const * filter_source  = object_type == HY_BL_DATASET_FILTER ? (_DataSetFilter const *)source_object : nil;
        _DataSet       const * dataset_source = filter_source ? nil : (_DataSet const *)source_object;

        switch (parameters.lLength) {
            case 2UL: { // get site->pattern map
                if (filter_source) {
                    receptacle->SetValue (new _Matrix (filter_source->duplicateMap),false,true, NULL);
                } else {
                    receptacle->SetValue (new _Matrix (dataset_source->DuplicateMap()),false,true, NULL);
                }
            }
            break;

            case 3UL: { // data parameters, or sequence string
                
                _String argument;
                try {
                    argument = _ProcessALiteralArgument (*GetIthParameter(2),current_program);
                } catch (const _String& err) {
                    //printf ("%s\n", err.get_str());
                    // not a string
                }
                if (argument.nonempty()) {
                    if (argument == kCharacters) {
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
                        receptacle->SetValue (new _Matrix (characters), false, true, NULL);
                    } else if (argument == kParameters) {
                        if (filter_source) {
                            _AssociativeList * parameterInfo = new _AssociativeList;

                            (*parameterInfo) < (_associative_list_key_value){"ATOM_SIZE", new _Constant (filter_source->GetUnitLength())}
                            < (_associative_list_key_value){"EXCLUSIONS", new _FString  (filter_source->GetExclusions())}
                            < (_associative_list_key_value){"SITES_STRING", new _FString  ((_String*)filter_source->theOriginalOrder.ListToPartitionString())}
                            < (_associative_list_key_value){"SEQUENCES_STRING", new _FString  ((_String*)filter_source->theNodeMap.ListToPartitionString())};

                            receptacle->SetValue (parameterInfo,false,true, NULL);

                        } else {
                            throw (argument.Enquote('\'') & " is only available for DataSetFilter objects");
                        }
                    } else if (argument == kConsensus) { // argument == _String("PARAMETERS")
                        if (filter_source) {
                            receptacle->SetValue (new _FString (new _String(filter_source->GenerateConsensusString())), false,true, NULL);
                        } else {
                            _DataSetFilter temp;
                            _SimpleList l1, l2;
                            temp.SetFilter (dataset_source, 1, l1, l2, false);
                            receptacle->SetValue (new _FString (new _String(temp.GenerateConsensusString())), false,true, NULL);
                        }
                    }
                } else {
                    long seqID = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix);

                    if (filter_source) {
                        if (seqID>=0 && seqID < filter_source->NumberSpecies()) {
                            receptacle->SetValue (new _FString (filter_source->GetSequenceCharacters(seqID)),false,true, NULL);
                        } else  if (seqID >= -4 && seqID <= -1) {
                            _SimpleList indices, map, counts;
                            _hy_dataset_filter_unique_match match_mode = kUniqueMatchExact;
                            switch (seqID) {
                                case -2:
                                    match_mode = kUniqueMatchExactOrGap;
                                    break;
                                case -3:
                                    match_mode = kUniqueMatchSuperset;
                                    break;
                                case -4:
                                    match_mode = kUniqueMatchPartialMatch;
                                    break;
                            }
                            long uniqueSequences = filter_source->FindUniqueSequences(indices, map, counts, match_mode);
                            _AssociativeList * parameterInfo = new _AssociativeList;
                            parameterInfo->MStore ("UNIQUE_SEQUENCES", new _Constant (uniqueSequences), false);
                            parameterInfo->MStore ("UNIQUE_INDICES",   new _Matrix   (indices), false);
                            parameterInfo->MStore ("SEQUENCE_MAP",     new _Matrix   (map), false);
                            parameterInfo->MStore ("UNIQUE_COUNTS",    new _Matrix   (counts), false);
                            receptacle->SetValue (parameterInfo,false,true, NULL);
                        }
                    } else { // filter_source
                        if (seqID>=0 && seqID < dataset_source->NoOfSpecies()) {
                            receptacle->SetValue (new _FString (dataset_source->GetSequenceCharacters(seqID)),false,true, NULL);
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
                                receptacle->SetValue (new _Constant (theValue),false,true, NULL);
                                DeleteObject     (res);
                            } else {
                                receptacle->SetValue (res,false,true, NULL);
                            }
                        } else {
                            bool count_gaps = hy_env::EnvVariableTrue(hy_env::harvest_frequencies_gap_options);
                            long filter_dimension = filter_source->GetDimension (true);

                            _Matrix * accumulator = new _Matrix (filter_dimension, 1, false, true),
                                    * storage     = new _Matrix (filter_dimension, 1, false, true);

                            _String *buffer = filter_source->MakeSiteBuffer();

                            for (long species_index = filter_source->NumberSpecies()-1; species_index >= 0; species_index --) {
                                filter_source->RetrieveState(site,species_index,*buffer, false);
                                filter_source->Translate2Frequencies (*buffer, storage->theData,  count_gaps);
                                *accumulator += *storage;
                            }
                            receptacle -> SetValue (accumulator, false,true, NULL);
                            BatchDelete(storage, buffer);

                        }
                    } else {
                        throw (_String ("Site index ") & site & " is invalid: must be in range " & "[0, " & (long)filter_source->GetPatternCount() & "]");
                    }
                } else {
                    throw _String("This set of arguments is only supported for DataSetFilter objects");
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

                        if (kPairwiseCountAmbiguitiesAverage == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingAverageFrequencyAware);
                        } else if (kPairwiseCountAmbiguitiesResolve == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolve);
                        } else if (kPairwiseCountAmbiguitiesSkip == *res_flag) {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingSkip);
                        } else {
                            res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolveFrequencyAware);
                        }

                        receptacle->SetValue (res,false,true, NULL);
                    } else {
                        throw (_String (seq1).Enquote() & "," & _String (seq2).Enquote() & " is an invalid sequence pair specification.");
                    }
                } else {
                    throw _String("This set of options is not supported for DataSet arguments");
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
                    
                    _List       modelPNames;
                 
                    PopulateAndSort ([&] (_AVLList & parameter_list) -> void  {
                         if (IsModelOfExplicitForm (object_index)) {
                             ((_Formula const*)source_object)->ScanFForVariables(parameter_list,false);
                         } else {
                             ((_Variable const*)source_object)->ScanForVariables(parameter_list,false);
                         }
                    }).Each ([&] (long value, unsigned long) -> void {
                        modelPNames << LocateVar(value)->GetName();
                    });
                   
                    result = new _Matrix (modelPNames);
                }
                break;
            }
        } else {
            _Variable* source_object = FetchVar(LocateVarByName (source_name));

            if (source_object && source_object->ObjectClass()==STRING) {
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
        receptacle->SetValue(result, false,true, NULL);
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleConstructCategoryMatrix (_ExecutionList& current_program) {

    static      _Trie        kRunOptions (
                                            _String ("COMPLETE"),_hyphyLFConstructCategoryMatrixConditionals,
                                            _String ("WEIGHTS"),_hyphyLFConstructCategoryMatrixWeights,
                                            _String ("SITE_LOG_LIKELIHOODS"), _hyphyLFConstructCategoryMatrixSiteProbabilities,
                                            _String ("CLASSES"), _hyphyLFConstructCategoryMatrixClasses,
                                            _String ("SHORT"), _hyphyLFConstructCategoryMatrixClasses
                                          );


    _Variable * receptacle = nil;
    current_program.advance();

    try {

        _Matrix*   result     = nil;
        receptacle = _ValidateStorageVariable (current_program);
        const _String source_name = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        long            object_type = HY_BL_LIKELIHOOD_FUNCTION  | HY_BL_TREE,
                        object_index;
        BaseRefConst    source_object = _GetHBLObjectByType (source_name, object_type, &object_index, &current_program);

        switch (object_type) {
            case HY_BL_LIKELIHOOD_FUNCTION: {
                _Matrix * partition_list  = nil;
                if (parameters.countitems () > 3) { // have a restricting partition
                    partition_list = (_Matrix*)_ProcessAnArgumentByType(*GetIthParameter(3), MATRIX, current_program, nil);
                }
                _SimpleList   included_partitions;
                _LikelihoodFunction * like_func = (_LikelihoodFunction * )source_object;

                like_func->ProcessPartitionList(included_partitions, partition_list);
                DeleteObject (partition_list);

                long run_mode = _hyphyLFConstructCategoryMatrixConditionals;

                if (parameters.countitems () > 2) {
                    long run_mode_long = kRunOptions.GetValueFromString(*GetIthParameter(2));
                    if (run_mode_long != kNotFound) {
                        run_mode = run_mode_long;
                    }
                }

                receptacle->SetValue(like_func->ConstructCategoryMatrix(included_partitions, run_mode ,true, receptacle->GetName()), false,true, NULL);
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
                        false,true, NULL
                    );
                }
            }
            break;
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
  
    _List   dynamic_variable_cleanup;

    try {
        _Matrix   * result     = nil;
        receptacle = _ValidateStorageVariable (current_program);
        _Matrix   * input_seqs = (_Matrix   *)_ProcessAnArgumentByType(*GetIthParameter(1), MATRIX, current_program, &dynamic_variable_cleanup);
 
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


        _AssociativeList* alignment_options  = (_AssociativeList   *)_ProcessAnArgumentByType(*GetIthParameter(2), ASSOCIATIVE_LIST, current_program, &dynamic_variable_cleanup);
        _FString        * char_vector        = (_FString*)          _EnsurePresenceOfKey (alignment_options, kCharacterMap, STRING);

        unsigned          long     char_count = 0UL;
        long              character_map_to_integers [256];
        InitializeArray(character_map_to_integers, 256, -1L);
        
        for (unsigned long cc = 0UL; cc < char_vector->get_str().length(); cc++) {
            unsigned char this_char = char_vector->get_str().get_uchar(cc);
            if (character_map_to_integers [this_char]>=0) {
                throw (_String ("Duplicate character ") & _String ((char)this_char).Enquote('\'') & " in " & kCharacterMap);
            } else {
                character_map_to_integers [this_char] = cc;
                char_count ++;
            }
        }
        if (char_count == 0) {
            throw _String("Null alphabet supplied");
        }

        bool        do_local       = _NumericValueFromKey (alignment_options,kGapLocal,0.0) > 0.5,
                    do_affine      = _NumericValueFromKey (alignment_options,kAffineGaps,0.0) > 0.5,
                    do_linear      = _NumericValueFromKey (alignment_options,kLinearSpace,1.0) > 0.5,
                    do_codon       = _NumericValueFromKey (alignment_options,kCodonAlign,0.0) > 0.5,
                    do_full_local  = do_codon && _NumericValueFromKey (alignment_options,kLocalAlignment,0.0) > 0.5;


        long        codon_count        = char_count * char_count * char_count,
                    expected_dimension = (do_codon ? codon_count : char_count) + 1UL;

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
            if (gap_c->get_str().length () != 1UL) {
                throw (_String ("Invalid gap character specification ") &  gap_c->get_str());
            }
            gap_character = gap_c->get_str().char_at(0UL);
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
        _String const * reference_sequence = & ((_FString*)input_seqs->GetFormula(0,0)->Compute())->get_str();


        for (unsigned long index2 = 1UL; index2 < input_seq_count; index2++) {
            _String const * sequence2 = & ((_FString*)input_seqs->GetFormula(0,index2)->Compute())->get_str();
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
                ops.list_data[reference_sequence->length()+1] = sequence2->length();
                ops.list_data[0]               = -1;

                score = LinearSpaceAlign(reference_sequence,sequence2,character_map_to_integers,score_matrix,
                                         gap_open,gap_extend,gap_open2,gap_extend2,
                                         do_local,do_affine,ops,score,0,
                                         reference_sequence->length(),0,sequence2->length(),buffers,0,alignment_route);

                delete[]    alignment_route;

                _StringBuffer     *result1 = new _StringBuffer (reference_sequence->length() + 1UL),
                                  *result2 = new _StringBuffer (size_allocation);

                long             last_column     = ops.list_data[ops.lLength-1];

                for (long position = (long)reference_sequence->length() - 1L; position>=0; position--) {
                    long current_column     = ops.list_data[position+1];

                    if (current_column<0) {
                        if (current_column == -2 /*|| (current_column == -3 && last_column == string2->sLength)*/) {
                            current_column = last_column;
                        } else if (current_column == -3) {
                            // find the next matched char or a -1
                            long    p   = position, s2p;
                            while ((ops.list_data[p+1]) < -1) {
                                p--;
                            }

                            s2p = ops.list_data[p+1];

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
                pairwise_alignment->MStore ("0", new _Constant (score), false);
                aligned_strings->MStore (_String((long)index2-1L), pairwise_alignment, false);
           } else { // not linear
                char * str1r = nil,
                     * str2r = nil;

               
                score = AlignStrings (reference_sequence->get_str() ? reference_sequence->get_str() : "",
                                      sequence2->get_str() ? sequence2->get_str() : "",
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
                    throw _String( "Internal Error in AlignStrings" );
                }
                pairwise_alignment->MStore ("0", new _Constant (score), false);
                aligned_strings->MStore (_String((long)index2-1L), pairwise_alignment, false);
            }
        }
        receptacle->SetValue(aligned_strings, false,true, NULL);
        
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
        BaseRefConst    source_object = _GetHBLObjectByType(*GetIthParameter(1), object_type, nil, &current_program);

        long      unit      = _ProcessNumericArgumentWithExceptions(*GetIthParameter(2),current_program.nameSpacePrefix),
                  atom      = _ProcessNumericArgumentWithExceptions(*GetIthParameter(3),current_program.nameSpacePrefix);

        bool      position_specific = _ProcessNumericArgumentWithExceptions(*GetIthParameter(4),current_program.nameSpacePrefix) > 0.5,
                  include_gaps       = hy_env::EnvVariableTrue(hy_env::harvest_frequencies_gap_options);

        switch (object_type) {
            case HY_BL_DATASET: {
                _String vertical_partition      (parameters.countitems () > 5 ? *GetIthParameter(5) : kEmptyString),
                        horizontal_partition    (parameters.countitems () > 6 ? *GetIthParameter(6) : kEmptyString);

                _DataSet const * dataset = (_DataSet const*)source_object;
                _SimpleList     processed_sequence_partition, processed_site_partition;
                dataset->ProcessPartition (horizontal_partition,processed_sequence_partition,false, 1);
                dataset->ProcessPartition (vertical_partition,processed_site_partition,true, 1);

                receptacle->SetValue (dataset->HarvestFrequencies(unit,atom,position_specific,processed_sequence_partition, processed_site_partition,include_gaps), false,true, NULL);
            }
            break;
            case HY_BL_DATASET_FILTER: {
                receptacle->SetValue (((_DataSetFilter const*)source_object)->HarvestFrequencies(unit,atom,position_specific,include_gaps), false,true, NULL);
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

        long       object_type = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM;
        _String    optimize_me = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);

        _LikelihoodFunction*   source_object = (_LikelihoodFunction*)_HYRetrieveBLObjectByNameMutable (AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix), object_type,nil,do_optimize==false);

        if (!source_object) { // Custom function (expression based)
            source_object = new _CustomFunction (optimize_me, current_program.nameSpacePrefix);
        }

        if (do_optimize) {
            if (parameters.countitems () > 2) { // have a restricting partition
                _List ref;
                receptacle -> SetValue(source_object->Optimize((_AssociativeList*)_ProcessAnArgumentByType(*GetIthParameter(2L), ASSOCIATIVE_LIST, current_program, &ref)),false,true, NULL);
            } else {
                receptacle -> SetValue(source_object->Optimize(),false,true, NULL);
            }
        } else {
            HBLObjectRef     covariance_parameters = hy_env::EnvVariableGet(hy_env::covariance_parameter, ASSOCIATIVE_LIST|STRING);
            _SimpleList   *restrictor = nil;
            switch (object_type) {
                case HY_BL_LIKELIHOOD_FUNCTION:
                case HY_BL_SCFG: {
                    if (covariance_parameters) { // only consider some variables
                        _SimpleList variable_ids;
                        if (covariance_parameters->ObjectClass () == ASSOCIATIVE_LIST) {
                            // a list of variables stored as keys in an associative array
                            _List*  restricted_variables = ((_AssociativeList*)covariance_parameters)->GetKeys();
                            for (unsigned long iid = 0; iid < restricted_variables->lLength; iid++) {
                                variable_ids << LocateVarByName (current_program.AddNameSpaceToID(*(_String*)(*restricted_variables)(iid)));
                            }
                            DeleteObject (restricted_variables);
                        } else { // STRING
                            variable_ids << LocateVarByName (current_program.AddNameSpaceToID(((_FString*)covariance_parameters)->get_str()));
                        }
                        if (!variable_ids.empty()) {
                            restrictor = new _SimpleList ();

                            for (unsigned long var_index = 0UL; var_index < variable_ids.lLength; var_index++) {
                                // TODO SLKP 20170706: this is a very inefficient linear search over and over again (in a large array)
                                long vID = source_object->GetIndependentVars().Find(variable_ids (var_index));
                                if (vID >= 0) (*restrictor) << vID;
                            }

                            if (restrictor->empty()) {
                                DeleteAndZeroObject (restrictor);
                            }
                        }
                    }
                }
                break;
                case HY_BL_BGM: {
                    _Matrix * bgm_cov = (_Matrix*)source_object->CovarianceMatrix(nil);
                    if (bgm_cov) {
                        receptacle->SetValue(bgm_cov,false,true, NULL);
                        return true;
                    } // TODO SLKP 20170706: handle the case when null is returned (why would that happen?); warn the user.
                }
                break;
            }

            receptacle->SetValue(source_object->CovarianceMatrix(restrictor),false,true, NULL);
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

    static const _String kLastSetOfConstraints    ("LAST_SET_OF_CONSTRAINTS");
    
    current_program.advance();
    
    _TreeIterator ** traversers = nil;
    unsigned long template_parameter_count = parameter_count() - 1L;
    deferSetFormula = new _SimpleList;

    auto cleanup = [template_parameter_count] (_TreeIterator ** t) -> void {
        if (t) {
            for (long k = 0L; k < template_parameter_count; k++) {
                if (t[k]) {
                    delete t[k];
                }
            }
            delete [] t;
        }
        delete deferSetFormula;
    };
    
    auto resolve_name = [] (_String const& my_name, BaseRefConst parent_name) -> _String const  {
        return my_name.Cut (((_String*)parent_name)->length() + 1, kStringEnd);
    };
    
    
    try {
        
        /** check the plug-in first **/
        
        _List templated_operations,
              parent_object_names;
        
        traversers = new _TreeIterator* [parameter_count()-1] {NULL};

        for (long k = 1L; k < parameter_count(); k++) {
            _CalcNode * this_object = (_CalcNode *)_CheckForExistingVariableByType (*GetIthParameter(k), current_program, TREE | TREE_NODE);
            if (this_object->ObjectClass () == TREE) {
                traversers[k-1] = new _TreeIterator ((_TheTree*)this_object, _HY_TREE_TRAVERSAL_POSTORDER);
                parent_object_names << this_object->GetName();
            } else {
                node<long> * cn = this_object->LocateMeInTree();
                if (!cn) {
                    throw (*GetIthParameter(k) & " is not a part of a tree object");
                }
                parent_object_names << this_object->ParentTree()->GetName();
                traversers[k-1] = new _TreeIterator (this_object, cn, _HY_TREE_TRAVERSAL_POSTORDER);
            }
            templated_operations < new _List;
        }
        
        /*
        _CalcNode * check = traversers[1]->Next();
        while (check) {
            printf ("%s\n", check->GetName()->get_str());
            check = traversers[1]->Next();
        }
        
        throw (_String("BAIL"));
        */

        _String  const  constraint_pattern = _ProcessALiteralArgument (*GetIthParameter(0), current_program);
        
        /** parse the contraint pattern using the standard exprssion parser
         because of '?' in template names, and things like __ which need to be deferred,
         the standard parser will fail with standard argument calls
         */

        _String parse_error, patttern_copy (constraint_pattern);
        _FormulaParsingContext fpc (&parse_error, current_program.nameSpacePrefix);
        fpc.allowTemplate() = '?';
        
        _Formula               pattern,
                               lhs;
        
        long    return_value = Parse (&pattern, patttern_copy, fpc, &lhs);
        
        if (return_value != HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT && return_value != HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT && return_value != HY_FORMULA_EXPRESSION) {
            throw (_String ("The template for ReplicateConstraint must be an assignment, with an expression based left-hand side, or simple expression"));
        }
        

        /** for 0 -- (K-1), where thisK is the largest "this" argument,
         stores the operation objects that reference these variables
         **/

        _Operation * lhs_op = nil;
        
        if (return_value != HY_FORMULA_EXPRESSION) {
            lhs_op = new _Operation (*LocateVar(fpc.assignmentRefID()));
            lhs.PushTerm (lhs_op);
            lhs_op->RemoveAReference();
        }
        
        _Formula * formulas [2] = {&lhs, &pattern};
        
        _List      substitution_variables (new _List, new _List);
        _List      list_of_constraints;
        
        auto       add_a_constraint = [return_value, &list_of_constraints] (_String* lhs, _String *rhs) -> void {
            _StringBuffer * new_constraint = new _StringBuffer;
            switch (return_value) {
                case HY_FORMULA_EXPRESSION:
                    new_constraint->AppendNewInstance(rhs);
                    break;
                case HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT:
                    (new_constraint->AppendNewInstance(lhs) << ":=").AppendNewInstance (rhs);
                    break;
                case HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT:
                    (new_constraint->AppendNewInstance(lhs) << "=").AppendNewInstance (rhs);
                    break;

            }
            
            list_of_constraints < new_constraint;
        };
        
        unsigned long reference_argument = 0UL; // either the first argument
                                                // or the argument that appears on the LHS
        
        _SimpleList   substitution_variable_by_index;
                        // this keeps track of which template parameter the i-th thisN.?.. argument comes from
        
        for (_Formula * f : formulas) {
            bool is_lhs = f == &lhs;
            f->GetList().ForEach([&templated_operations, is_lhs, &reference_argument, &substitution_variables, &substitution_variable_by_index] (BaseRef op, unsigned long index) -> void {
                _Variable const * operation_var = ((_Operation*) op)->RetrieveVar();
                if (operation_var) {
                    _SimpleList pattern_match (operation_var->GetName()->RegExpMatch(hy_replicate_constraint_regexp, 0));
                    if (pattern_match.nonempty()) {
                        unsigned long var_index = operation_var->GetName()->Cut (pattern_match (2), pattern_match (3)).to_long() - 1UL;
                        if (var_index >= templated_operations.countitems()) {
                            throw (operation_var->GetName()->Enquote() & " does not have a matched positional argument");
                        }
                        if (is_lhs) {
                            reference_argument = var_index;
                        }
                        _List * term_record = new _List;
                        term_record->AppendNewInstance (new _String (*operation_var->GetName(), pattern_match (4), pattern_match (5)));
                        (*term_record) << op;
                        //printf ("[term_record] %s\n", ((_String*) (term_record->toStr()))->get_str());

                        ((_List*)templated_operations.GetItem(var_index))->AppendNewInstance(term_record);
                        *((_List*)(substitution_variables.GetItem(0))) << operation_var->GetName();
                        substitution_variable_by_index << var_index;
                    }
                }
            });
        }
        
        
        templated_operations.ForEach ([this] (BaseRef items, long index) -> void {
            if (((_List*)items)->empty()) {
                throw (GetIthParameter(index + 1L)->Enquote() & " (argument " & _String (index+1L) & ") did not appear in the contstraint expression");
            }
        });
        
        /** perform a joint traversal of all the template objects
            the first template object will be used to create the list
            of local parameters that could be used to populate template contraints
         */
        
        
        auto incompatible_error = [this, reference_argument] (unsigned long i, const _String& name) -> void {
            throw (GetIthParameter(i + 1L)->Enquote() & " (argument " & _String ((long)i+1L) &
                   ") is topologically incompatible with the reference argument " & GetIthParameter(reference_argument + 1UL)->Enquote()) & " at node " & name.Enquote();
        };
        
        for (;;) {
            _CalcNode * reference_iteratee = traversers[reference_argument]->Next();
            node <long> * reference_node   = traversers[reference_argument]->GetNode();
            
            
            if (reference_iteratee) {
                //printf ("%s\n", reference_iteratee->GetName()->get_str());
                _List iteratees;
                iteratees << reference_iteratee;
                for (unsigned long i = 0UL; i < template_parameter_count; i++) {
                    if (reference_argument != i) {
                        _CalcNode * current_iteratee = traversers[i]->Next();
                        if (!current_iteratee) {
                            incompatible_error (i, *reference_iteratee->GetName());
                        }
                        //printf ("%s\n", current_iteratee->GetName()->get_str());
                        node <long> *current_node = traversers[i]->GetNode();
                        if (current_node->get_num_nodes() != reference_node->get_num_nodes()) {
                            incompatible_error (i,*reference_iteratee->GetName());
                        }
                        iteratees << current_iteratee;
                    }
                }
                
                auto handle_variable = [&parent_object_names,resolve_name] (_Variable* local_k, _List * conditions, _List * parameter_set, _List *matched_subexpressions, long condition_index) -> void {
                    _String local_name = resolve_name (*local_k->GetName(), parent_object_names.GetItem(condition_index));
                    _SimpleList matched_conditions;
                    
                    if (conditions->Every ([&local_name,&matched_conditions] (long condition, unsigned long i) -> bool {
                        _String const * match_pattern = (_String*)((_List*)condition)->GetItem (0);
                        //printf ("[name] %s [pattern] %s\n",local_name.get_str(), match_pattern->get_str());
                        if (i == 0) {
                            return local_name.EqualWithWildChar (*match_pattern,'?', 0, 0, &matched_conditions);
                            
                        } else {
                            return local_name.EqualWithWildChar (*match_pattern,'?');
                        }
                    })) {
                        //printf ("[%ld] %s (%s, %s)\n", condition_index,local_k->GetName()->get_str(), local_name.get_str(), ((_String*)parent_object_names.GetItem(condition_index))->get_str());
                        *parameter_set << local_k;
                        _List * subexpressions = new _List;
                        for (long k = 0L; k < matched_conditions.countitems(); k+=2) {
                            subexpressions->AppendNewInstance (new _String (local_name, matched_conditions(k), matched_conditions (k+1)));
                        }
                        //ObjectToConsole(subexpressions);
                        *matched_subexpressions < subexpressions;
                    }
                };

                if (reference_iteratee->HasLocals()) { // stuff to do
                    _List parameter_sets,
                          matched_subexpressions;
                    for (unsigned long i = 0UL; i < template_parameter_count; i++) {
                        parameter_sets < new _List;
                        matched_subexpressions < new _List;
                    }
                    

                    /** now see if there any applicable contstraints that can be generated
                     if there is a LHS expression, then only independent variables will be pulled
                     */
                    
                    long independent_count = reference_iteratee->CountIndependents();
                    _List * conditions = ((_List*)templated_operations.GetItem(reference_argument));
                    _List * reference_parameters = (_List*)parameter_sets.GetItem(reference_argument);
                    _List* matched_reference_subexpressions = (_List*)matched_subexpressions.GetItem(reference_argument);
                        /*
                            this will store, for each reference argument, the parts of parameter names
                            that were matched by the '?' in the template expression.
                            The convention is that the first set takes precedence, so that
                            this1.?.? := this1.?.synRate will store the matches for the two '?' from the first expression
                         */
                    
                    
                    for  (long k = 0L; k < independent_count; k++) {
                        handle_variable ( reference_iteratee->GetIthIndependent(k), conditions, reference_parameters, matched_reference_subexpressions, reference_argument);
                    }
                    
                    if (!lhs_op) {
                        for  (long k = 0L; k < reference_iteratee->CountDependents(); k++) {
                            handle_variable ( reference_iteratee->GetIthDependent(k), conditions, reference_parameters, matched_reference_subexpressions, reference_argument);
                        }
                    }
                    
                    /** now check to see if subsequent template arguments fit reference template parameters
                        If there are multiple parameter matches for any subsequent template, then
                     */
                    
                    if (reference_parameters->nonempty()) {
                        for (unsigned long i = 0UL; i < template_parameter_count; i++) {
                            if (i != reference_argument) {
                                _CalcNode  * template_iteratee = (_CalcNode*)iteratees.GetItem(i);
                                _List      * template_variables = (_List*)parameter_sets.GetItem(i);
                                _List      * matched_template_subexpressions = (_List*)matched_subexpressions.GetItem(i);
                                _List      * template_conditions = ((_List*)templated_operations.GetItem(i));
                                long      variable_count = template_iteratee->CountIndependents();
                                
                                
                                // see which parameters match the template conditions
                                for  (long k = 0L; k < variable_count; k++) {
                                    handle_variable ( template_iteratee->GetIthIndependent(k), template_conditions, template_variables, matched_template_subexpressions, i);
                                }
                                variable_count = template_iteratee->CountDependents();
                                for  (long k = 0L; k < variable_count; k++) {
                                    handle_variable ( template_iteratee->GetIthDependent(k), template_conditions, template_variables, matched_template_subexpressions, i);
                                }
                            }
                        }
                        // now for each reference parameter, see if it is matched by the corresponding template paramters
                        
                        reference_parameters->ForEach([&substitution_variables, &pattern, &lhs, &matched_subexpressions, reference_argument, &parameter_sets, add_a_constraint, &current_program, &substitution_variable_by_index] (BaseRef  ref, unsigned long i) -> void {
                            _Variable * ref_var = (_Variable *)ref;
                            _List * reference_subexpressions    = (_List*)matched_subexpressions.GetItem(reference_argument,i);
                            _List   local_substitution_variables;
                            
                            //ObjectToConsole(reference_subexpressions);
                            
                            if   (parameter_sets.Every ([reference_argument, ref_var, &local_substitution_variables, matched_subexpressions, reference_subexpressions] (long template_i, unsigned long i2) -> bool {
                                if (i2 != reference_argument) {
                                    _List * template_parameter_set = (_List *)template_i;
                                    _List * matched_template_subexpressions = (_List*)matched_subexpressions.GetItem(i2);
                                    
                                    return template_parameter_set->Any ([matched_template_subexpressions, reference_subexpressions, &local_substitution_variables] (long argument_ij, unsigned long i3) -> bool {
                                        _List *  check_subs = (_List*) matched_template_subexpressions->GetItem (i3);
                                        for (long i = 0; i < Minimum(check_subs->countitems(), reference_subexpressions->countitems()); i++) {
                                            if (*(_String*)check_subs->GetItem(i) != *(_String*)reference_subexpressions->GetItem(i) ) {
                                               // printf ("[fail check] %s %s [%d]\n", ((_String*)check_subs->GetItem(i))->get_str(),((_String*)reference_subexpressions->GetItem(i))->get_str(), check_subs->countitems());
                                                return false;
                                            }
                                        }
                                        local_substitution_variables << (_Variable *)argument_ij;
                                        return true;
                                    });
                                } else {
                                    local_substitution_variables << ref_var;
                                    return true;
                                }
                            })) { // this template argument is OK to go
                                ((_List*)(substitution_variables.GetItem(1)))->Duplicate(substitution_variables.GetItem(0));
                                //printf ("%s\n", ((_String*)local_substitution_variables.toStr())->get_str());
                                substitution_variable_by_index.Each ([&substitution_variables,&local_substitution_variables] (long template_index, unsigned long idx) -> void {
                                    _String * sub_name = ((_Variable*) local_substitution_variables.GetItem(template_index))->GetName();
                                    ((_List*)(substitution_variables.GetItem(1)))->InsertElement (sub_name, idx, false, true);
                                });
                                

                                /*local_substitution_variables.ForEach([&substitution_variables] (BaseRef v, unsigned long) -> void {
                                     *((_List*)(substitution_variables.GetItem(1))) << ((_Variable*) v)->GetName();
                                });*/
                                
                                add_a_constraint (lhs.IsEmpty() ? nil : ((_String*)lhs.toStr (kFormulaStringConversionNormal, &substitution_variables)),
                                                  ((_String*)pattern.toStr (kFormulaStringConversionNormal, &substitution_variables)));
                                // apply the constraint
                            }
                            
                        });
                    }

                }
            } else { // possible termination
                for (unsigned long i = 0UL; i < template_parameter_count; i++) {
                    if (reference_argument != i) {
                        if (traversers[i]->Next()) {
                            incompatible_error (i, "Root node");
                        }
                    }
                }
                break;
            }
        }
        
        list_of_constraints.ForEach ([&current_program] (BaseRefConst constraint_i, unsigned long) -> void {
            _String constraint = *(_String*)constraint_i;
            _Formula rhs, lhs;
            _FormulaParsingContext fpc (nil, current_program.nameSpacePrefix);
            long result = Parse (&rhs,constraint,fpc,&lhs);
            ExecuteFormula(&rhs,&lhs,result,fpc.assignmentRefID(),current_program.nameSpacePrefix,fpc.assignmentRefType());
        });

        
        _StringBuffer * generated_constraints = new _StringBuffer (list_of_constraints.Join (";\n"));
        if (generated_constraints->nonempty()) {
            *generated_constraints << ";";
        }
        
        hy_env::EnvVariableSet(kLastSetOfConstraints, new _FString (generated_constraints), false);
        
        FinishDeferredSF ();
        
        // first find the pattern


    } catch (const _String& error) {
        cleanup (traversers);
        return  _DefaultExceptionHandler (nil, error, current_program);
    }
    
    cleanup (traversers);
    return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleComputeLFFunction (_ExecutionList& current_program) {

  const static _String kLFStartCompute ("LF_START_COMPUTE"),
                       kLFDoneCompute  ("LF_DONE_COMPUTE"),
                       kLFTrackCache   ("LF_TRACK_CACHE"),
                       kLFAbandonCache ("LF_ABANDON_CACHE");

  current_program.advance();
  _Variable * receptacle = nil;

  try {


    _String    const op_kind = * GetIthParameter(1UL);

    long       object_type = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM;
    _LikelihoodFunction*    source_object = (_LikelihoodFunction*)_GetHBLObjectByType(AppendContainerName (*GetIthParameter(0UL), current_program.nameSpacePrefix),object_type, nil,&current_program);

    if (op_kind == kLFStartCompute) {
      source_object->PrepareToCompute(true);
     } else if (op_kind == kLFDoneCompute) {
      source_object->FlushLocalUpdatePolicy();
      source_object->DoneComputing (true);
    } else {
      if (!source_object->HasBeenSetup()) {
        throw (_String("Please call LFCompute (, ") & *GetIthParameter (0UL)& kLFStartCompute & ") before evaluating the likelihood function");
      } else {
        if (op_kind == kLFTrackCache) {
          source_object->DetermineLocalUpdatePolicy();
        } else if (op_kind == kLFAbandonCache) {
          source_object->FlushLocalUpdatePolicy();
        } else {
          receptacle = _ValidateStorageVariable (current_program, 1UL);
          receptacle->SetValue (new _Constant (source_object->Compute()), false,true, NULL);
        }
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

    _Variable *    source_model = (_Variable*)_HYRetrieveBLObjectByNameMutable (source_name, object_type,&model_index,false);

    if (!source_model && raw_model_name != kUseNoModel) {
        throw (source_name.Enquote() & " does not refer to a valid defined substitution model and is not " & kUseNoModel);
    }
      
    lastMatrixDeclared = model_index;

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }

  return true;
}


//____________________________________________________________________________________

bool      _ElementaryCommand::HandleInitializeIterator (_ExecutionList& current_program) {
/*
    parameters:
        [0] -> the string ID of the object to iterate over
        [1] -> the object to iteratre over
 
    
    simpleParameters:
        [0] -> command index for the corresponding advance command
        [1] -> type of object; MATRIX, etc
        [2] -> if dict, then this stores the pointer to AVLListXLIterator
            -> if tree, then this stores the pointer to a tree iterator
*/
  current_program.advance();
  try {
      // clean up previous iterator states
      if (parameters.countitems() > 1) {
          parameters.Delete (1L);
          if (simpleParameters.get (1) == ASSOCIATIVE_LIST) {
              delete ((AVLListXLIterator*)simpleParameters.get(2));
          } else {
              if (simpleParameters.get (1) == TREE || simpleParameters.get (1) == TOPOLOGY) {
                  delete ((node_iterator<long>*)simpleParameters.get(2));
              }
          }
       }
      
      if (simpleParameters.countitems() > 2) {
          simpleParameters.Delete (2);
      }
      
      
      HBLObjectRef iterator_substrate = _ProcessAnArgumentByType (*GetIthParameter(0UL), ASSOCIATIVE_LIST | MATRIX | TOPOLOGY | TREE, current_program, &parameters);
      
      _ElementaryCommand * advance_command = current_program.GetIthCommand(simpleParameters.get (0));
      
      if (!advance_command || advance_command->code != HY_HBL_COMMAND_ADVANCE_ITERATOR) {
          throw _String ("Iterator init command is not linked with an advance command");
      }
      
      advance_command->parameters.DeleteTail(advance_command->simpleParameters.get (1));
      
      if (iterator_substrate->ObjectClass() == ASSOCIATIVE_LIST) {
          _AssociativeList* source_object = (_AssociativeList*)iterator_substrate;
          simpleParameters[1] = ASSOCIATIVE_LIST; // indicate that this is an ASSOCIATIVE_LIST
          simpleParameters << (long) new AVLListXLIterator (source_object->ListIterator ());
      } else { // MATRIX
          if (iterator_substrate->ObjectClass() == MATRIX) {
               simpleParameters[1] = MATRIX; // indicate that tjis is a MATRIX
          } else {
              simpleParameters[1] = iterator_substrate->ObjectClass();
              _TreeTopology * source_tree = (_TreeTopology*)iterator_substrate;
              simpleParameters << (long) new node_iterator<long> (&source_tree->GetRoot(), _HY_TREE_TRAVERSAL_POSTORDER);
          }
      }
  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }

  return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleAdvanceIterator(_ExecutionList& current_program) {
 /*
     parameters:
         [0-2] -> Id of the variable to store stuff in
         [+1] -> the source object
         [+2] -> _Variable to store the value in
         [+3] -> _Variable to store key/index in (could be nil)
         [+4] -> _Variable to store column index in (if the object is _Matrix)
         
  
     
     simpleParameters:
         [0] -> index of the iterator init command in current_program
         [1] -> number of receptable variables
         [2] -> for _Matrix iterator, the row index
         [3] -> for _Matrix iterator, the column index
*/
  current_program.advance();
  try {
      
      long reciever_count = simpleParameters.get (1);
      
      _ElementaryCommand * init_command = current_program.GetIthCommand(simpleParameters.get (0));
      if (!init_command || init_command->code != HY_HBL_COMMAND_INIT_ITERATOR) {
          throw (_String ("Iterator advance command is not linked with an iterator initiaizer"));
      }
      
      long source_object_class = init_command->simpleParameters.get (1);
      
      bool first_in = false;
      if (source_object_class != MATRIX && reciever_count == 3) {
          throw (_String ("Iterators over dictionaries do not support the 3 argument form"));
      }
      
      if (parameters.countitems () == reciever_count || parameters.GetItem(reciever_count) != init_command->parameters.GetItem(1)) {
          parameters.DeleteTail (reciever_count, true);
          parameters << init_command->parameters.GetItem(1);
          first_in = true;
          for (long k = 0L; k < reciever_count; k++) {
              parameters << _ValidateStorageVariable (current_program, k);
          }
          if (source_object_class == MATRIX) {
              if (simpleParameters.countitems() == 2) {
                  simpleParameters << 0 << 0;
              } else {
                  simpleParameters[2] = 0;
                  simpleParameters[3] = 0;
              }
          }
      }
      
      if (source_object_class == MATRIX) {
          _Matrix* source_object = (_Matrix*)parameters.GetItem (reciever_count);
          long row    = simpleParameters.get (2),
               column = simpleParameters.get (3);
          
          if (!first_in) {
              column ++;
              if (column >= source_object->GetVDim()) {
                  column = 0;
                  row ++;
              }
           }
          
           if (row >= source_object->GetHDim()) {
               ((_Variable*)parameters.GetItem (reciever_count << 1))->SetValue (new _FString (kEndIteration), false, false, NULL); // end loop
           } else {
               ((_Variable*)parameters.GetItem (reciever_count << 1))->SetValue (source_object->GetMatrixCell(row, column), false, false, NULL);
               if (reciever_count == 2) {
                   ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue (new _Constant (row*source_object->GetVDim () + column), false, false, NULL);
               } else if (reciever_count == 3) {
                   ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue (new _Constant (row), false, false,NULL);
                   ((_Variable*)parameters.GetItem (reciever_count+2))->SetValue (new _Constant (column ), false, false,NULL);
               }
           }
            
          simpleParameters[2] = row;
          simpleParameters[3] = column;


          
      } else {
          if (source_object_class == ASSOCIATIVE_LIST) {
              AVLListXLIterator * it = (AVLListXLIterator *)init_command->simpleParameters.get (2);
              if (first_in) {
                  it->begin();
              } else {
                  ++(*it);
              }
              if (!it->is_done()) { // iterator not finished
                  AVLListXLIteratorKeyValue state = *(*it);
                  state.get_object()->AddAReference();
                  if (reciever_count > 1) {
                      ((_Variable*)parameters.GetItem (reciever_count+2))->SetValue ((HBLObjectRef)state.get_object(), false, false, NULL);
                      ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue (new _FString (*state.get_key()), false, false,NULL);
                  } else {
                      ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue ((HBLObjectRef)state.get_object(), false, false,NULL);
                  }
              } else {
                  ((_Variable*)parameters.GetItem (reciever_count << 1))->SetValue (new _FString (kEndIteration), false, false,NULL);
                  // iterator done
              }
          } else {
              node_iterator<long> * tree_iterator = (node_iterator<long> *)init_command->simpleParameters.get (2);
              if (first_in) {
                  if (simpleParameters.countitems() < 3) {
                      simpleParameters << 0;
                  } else {
                      simpleParameters[2] = 0;
                  }
              }
              node<long>*topTraverser = tree_iterator->Next();
              if (topTraverser && !topTraverser->is_root()) {
                  _FString *node_name;
                  if (source_object_class == TREE)
                      node_name = new _FString (map_node_to_calcnode (topTraverser)->ContextFreeName());
                  else
                      node_name = new _FString (((_TreeTopology*)parameters.GetItem (1))->GetNodeName(topTraverser));
                        
                  if (reciever_count > 1) {
                      ((_Variable*)parameters.GetItem (reciever_count+2))->SetValue (node_name, false, false,NULL);
                      ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue (new _Constant (simpleParameters.get(2)), false, false,NULL);
                  } else {
                      ((_Variable*)parameters.GetItem (reciever_count+1))->SetValue (node_name, false, false,NULL);
                  }
                  simpleParameters[2]++;
              } else {
                ((_Variable*)parameters.GetItem (reciever_count << 1))->SetValue (new _FString (kEndIteration), false, false,NULL);
              }
              
          }
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
      
          auto throw_error = [&] (void) -> void {
              throw _String ("Current script requires at least version ")& requested_version &" of HyPhy. Had version " &kHyPhyVersion & ". Please download an updated version from http://www.hyphy.org or github.com/veg/hyphy and try again.";
          };
      
        const _List   local_version    = kHyPhyVersion.Tokenize ("."),
                      required_version = requested_version.Tokenize(".");
     
          unsigned long const  upper_bound = MIN (local_version.countitems (), required_version.countitems());
      
      
      for (unsigned long i = 0UL; i < upper_bound; i++) {
          hyFloat local_number = ((_String*)local_version.GetItem(i))->to_float(),
                  required_number = ((_String*)required_version.GetItem(i))->to_float();
          
          if (local_number > required_number) {
              return true;
          }
          if (local_number < required_number) {
              throw_error();
          }
      }

      if (required_version.countitems() <= upper_bound) {
          return true;
      }
      throw_error ();

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
}


    
  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleDeleteObject(_ExecutionList& current_program){
    
  current_program.advance();

  const static _String kShallow (":shallow");
    
  bool do_shallow = parameter_count() >= 2 || *GetIthParameter(parameter_count()-1, false) == kShallow;

  for (unsigned long i = 0UL; i < parameter_count(); i++) {
    long       requested_type = HY_BL_LIKELIHOOD_FUNCTION,
               object_index   = kNotFound;
    BaseRef    source_object = _HYRetrieveBLObjectByNameMutable (AppendContainerName(*GetIthParameter(i),current_program.nameSpacePrefix), requested_type,&object_index,false);

    if  (source_object) {
      KillLFRecord (object_index,!do_shallow);
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
        receptacle->SetValue(new _FString (url,false),false,true,NULL);
      } else {
        throw (_String ("Could not fetch ") & url.Enquote());
      }
    } else {
      if (*action == save_to_file_action) {
        _String file_name = _ProcessALiteralArgument (*GetIthParameter(1UL),current_program);
        if (!ProcessFileName (file_name, true,true,(hyPointer)current_program.nameSpacePrefix,false,&current_program)) {
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
  static _String last_model_used, kPromptText ("Select a standard model");

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
      _DataSetFilter const *    source_filter = (_DataSetFilter const*)_GetHBLObjectByType(source_name, object_type, nil, &current_program);


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
      bool need_to_prompt_user = true;

      if (current_program.has_stdin_redirect() || current_program.has_keyword_arguments()) {
        _String option;
        bool    has_input = true;
        try {
            _FString * redirect = (_FString*)hy_env::EnvVariableGet(hy_env::fprintf_redirect, STRING);
            option = current_program.FetchFromStdinRedirect (&kPromptText, false, !(redirect&&redirect->has_data()));
        } catch (const _String& e) {
            if (e != kNoKWMatch) {
                throw (e);
            }
            has_input = false;
        } 
        if (has_input) {
            model_id = matching_models.FindOnCondition( [&] (long index, unsigned long) -> bool {
              return option == * (_String*) templateModelList.GetItem (index,0);
            });


            if (model_id == kNotFound) {
              throw (option.Enquote() & " is not a valid model (with input redirect)");
            }
            
            need_to_prompt_user = false;
        }
      }
          
      if (need_to_prompt_user) {
        #ifdef __HEADLESS__
          throw _String("Unhandled standard input interaction in SelectTemplateModel for headless HyPhy");
        #endif



        for (int i = 0; i < kMaxDialogPrompts; i++) {
          printf ("\n\n               +--------------------------+\n");
          printf (    "               | %s. |\n", kPromptText.get_str());
          printf (    "               +--------------------------+\n\n\n");

          for (model_id = 0; model_id<matching_models.lLength; model_id++) {
            printf ("\n\t(%s):%s",((_String*)templateModelList.GetItem(matching_models(model_id),0))->get_str(),
                                  ((_String*)templateModelList.GetItem(matching_models(model_id),1))->get_str());
          }
          printf ("\n\n Please type in the abbreviation for the model you want to use:");
          _String const user_choice = StringFromConsole();

          model_id = matching_models.FindOnCondition( [&] (long index, unsigned long) -> bool {
            return user_choice.EqualIgnoringCase(*(_String*) templateModelList.GetItem (index,0));
          });

          if (model_id != kNotFound) {
            break;
          }
        }

        if (model_id == kNotFound) {
          throw _String("Dialog did not return a valid choice after maximum allowed number of tries");
          return false;
        }

      }

      _String  model_file = GetStandardDirectory (HY_HBL_DIRECTORY_TEMPLATE_MODELS) & *(_String*)templateModelList.GetItem(matching_models(model_id),4);

      _ExecutionList       std_model;
      PushFilePath        (model_file, false);
      ReadBatchFile       (model_file,std_model);
      PopFilePath         ();
      last_model_used    = model_file;
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

bool      _ElementaryCommand::HandleMPIReceive (_ExecutionList& current_program){
  current_program.advance();
  _Variable * receptacle = nil;

  try {

#ifdef __HYPHYMPI__

    receptacle = _ValidateStorageVariable (current_program, 2UL);
    _Variable* node_index_storage = _ValidateStorageVariable (current_program, 1UL);

    long target_node = _ProcessNumericArgumentWithExceptions(*GetIthParameter(0UL), current_program.nameSpacePrefix),
    node_count  = hy_env::EnvVariableGetNumber(hy_env::mpi_node_count);

    if (target_node < -1L || target_node >= node_count) {
      throw (GetIthParameter(1UL)->Enquote () & " (=" & node_count & ") must be a valid MPI node index (or -1 to accept from any node");
    }

    long received_from;
    receptacle->SetValue(new _FString (MPIRecvString (target_node,received_from)), false, true,NULL);
    node_index_storage->SetValue (new _Constant (received_from), false, true, NULL);

#else
    throw _String("Command not supported for non-MPI versions of HyPhy. HBL scripts need to check for MPI before calling MPI features");
#endif

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (receptacle, error, current_program);
  }

  return true;
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleMPISend (_ExecutionList& current_program){
  current_program.advance();
  
  _List dynamic_variable_manager;
  
  try {

#ifdef __HYPHYMPI__
   long target_node = _ProcessNumericArgumentWithExceptions(*GetIthParameter(0UL), current_program.nameSpacePrefix),
         node_count  = hy_env::EnvVariableGetNumber(hy_env::mpi_node_count);

    if (target_node < 0L || target_node >= node_count) {
        throw (GetIthParameter(1UL)->Enquote () & " (=" & node_count & ") is not a valid MPI node index; valud range is " & target_node & " to " & (node_count-1));
    }

    _StringBuffer message_to_send (1024UL);

    if (parameter_count() > 2) { // this is the case of MPISend (node, filepath, {input option}
      _AssociativeList * arguments = (_AssociativeList*)_ProcessAnArgumentByType(*GetIthParameter(2), ASSOCIATIVE_LIST, current_program, &dynamic_variable_manager);
      _String file_path = *GetIthParameter(1UL);
      if (! ProcessFileName(file_path,false,true,(hyPointer)current_program.nameSpacePrefix)) {
        throw (GetIthParameter(1UL)->Enquote() & " is an ivalid file path");
      }
      (message_to_send << _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_EXECUTE_A_FILE) << file_path.Enquote() << ',')
        .AppendNewInstance(arguments->Serialize(0UL)) << ");";
    } else {
        // is this a likelihood function?
        long type = HY_BL_LIKELIHOOD_FUNCTION;
        try {
          ((_LikelihoodFunction*) _GetHBLObjectByTypeMutable(AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix), type))->SerializeLF(message_to_send, _hyphyLFSerializeModeOptimize);
        } catch (const _String &) {
            // catch literal cases here
            message_to_send = _ProcessALiteralArgument(*GetIthParameter(1UL), current_program);
        }
    }

    if (message_to_send.nonempty()) {
      MPISendString(message_to_send, target_node);
    } else {
      throw (_String ("An invalid (empty) MPI message"));
    }
#else
    throw _String("Command not supported for non-MPI versions of HyPhy. HBL scripts need to check for MPI before calling MPI features");
#endif

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }

  return true;
}



//____________________________________________________________________________________

bool      _ElementaryCommand::HandleSetParameter (_ExecutionList& current_program) {

  static const _String kBGMNodeOrder ("BGM_NODE_ORDER"),
                       kBGMGraph     ("BGM_GRAPH_MATRIX"),
                       kBGMScores    ("BGM_SCORE_CACHE"),
                       kBGMConstraintMx ("BGM_CONSTRAINT_MATRIX"),
                       kBGMParameters   ("BGM_NETWORK_PARAMETERS");

  current_program.advance();
  
  _List dynamic_variable_manager;

  try {
    _String object_to_change = *GetIthParameter(0UL);

    /* handle special cases */
    if (object_to_change == hy_env::random_seed) {
      hy_env::EnvVariableSet(hy_env::random_seed, new _Constant (hy_random_seed = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix) ), false);
      init_genrand(hy_random_seed);
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

    if (object_to_change == hy_env::execution_mode) {
      current_program.errorHandlingMode = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix);
      return true;
    }

    if (object_to_change == hy_env::status_bar_update_string) {
      SetStatusLineUser (_ProcessALiteralArgument (*GetIthParameter(1), current_program));
      return true;
    }
      

    const _String source_name   = AppendContainerName (*GetIthParameter(0), current_program.nameSpacePrefix);

    long          object_type = HY_BL_ANY,
                  object_index;

    BaseRef       source_object;
    _String set_this_attribute = *GetIthParameter(1UL);

    try {
        source_object = _GetHBLObjectByTypeMutable (source_name, object_type, &object_index);
    } catch (const _String& error) { // handle cases when the source is not an HBL object
        
        _CalcNode* tree_node = nil;
        

        if (source_name.IsValidIdentifier(fIDAllowFirstNumeric|fIDAllowCompound)) {
             tree_node = (_CalcNode*)FetchObjectFromVariableByType(&source_name, TREE_NODE);
        } else{
            _String converted = source_name.ConvertToAnIdent(fIDAllowFirstNumeric|fIDAllowCompound);
            tree_node = (_CalcNode*)FetchObjectFromVariableByType(&converted, TREE_NODE);
        }
        
      if (tree_node) {
        if (set_this_attribute == _String("MODEL")) {
          _String model_name = AppendContainerName(*GetIthParameter(2UL),current_program.nameSpacePrefix);
          long model_type = HY_BL_MODEL, model_index;
          _Matrix* model_object            = (_Matrix*)_GetHBLObjectByTypeMutable(model_name, model_type, &model_index);
          _TheTree * parent_tree = (_TheTree * )tree_node->ParentTree();
          if (!parent_tree) {
            throw (GetIthParameter(0UL)->Enquote() & " is an orphaned tree node (the parent tree has been deleted)");
          }
          tree_node->ReplaceModel (model_name, parent_tree);
          long partition_id, likelihood_function_id = ((_TheTree*)parent_tree->Compute())->IsLinkedToALF(partition_id);
          if (likelihood_function_id>=0){
            //throw(parent_tree->GetName()->Enquote() & " is linked to a likelihood function (" & *GetObjectNameByType (HY_BL_LIKELIHOOD_FUNCTION, likelihood_function_id) &") and cannot be modified ");
            //return false;
            ((_LikelihoodFunction*)likeFuncList (likelihood_function_id))->Rebuild(true);
          }
        } else {
          throw  (set_this_attribute.Enquote() & " is not a supported parameter type for a tree node argument");
        }
      } else {
        throw (GetIthParameter(0UL)->Enquote() & " is not a supported object type");
      }
      return true;
    }


    switch (object_type) {
      case HY_BL_BGM: { // BGM Branch

        _BayesianGraphicalModel * bgm = (_BayesianGraphicalModel *) source_object;
        long    num_nodes = bgm->GetNumNodes();


          // set data matrix
        if (set_this_attribute == kBGMData) {
          _Matrix     * data_mx = (_Matrix *) _ProcessAnArgumentByType(*GetIthParameter(2UL), MATRIX, current_program, &dynamic_variable_manager);

            if (data_mx->GetVDim() == num_nodes) {
              bgm->SetDataMatrix (data_mx);
            } else {
              throw (_String("Data matrix columns (") & data_mx->GetVDim() & " ) does not match number of nodes in graph (" & num_nodes & ")");
            }

        }
        else if (set_this_attribute == kBGMScores) {
          bgm->ImportCache((_AssociativeList *)_ProcessAnArgumentByType(*GetIthParameter(2UL), ASSOCIATIVE_LIST, current_program, &dynamic_variable_manager));
        } // set structure to user-specified adjacency matrix
        else if (set_this_attribute == kBGMGraph) {
          _Matrix     * graphMx   = (_Matrix *) _ProcessAnArgumentByType(*GetIthParameter(2UL), MATRIX, current_program, &dynamic_variable_manager);

          if (graphMx->check_dimension(num_nodes, num_nodes)) {
            bgm->SetStructure ((_Matrix *) graphMx->makeDynamic());
          } else {
            throw _String("Dimension of graph does not match current graph");
          }
        } // set constraint matrix
        else if (set_this_attribute == kBGMConstraintMx) {
          _Matrix     * constraint_mx  = (_Matrix *) _ProcessAnArgumentByType(*GetIthParameter(2UL), MATRIX, current_program, &dynamic_variable_manager);
          if (constraint_mx->check_dimension(num_nodes, num_nodes)) {
            bgm->SetConstraints ((_Matrix *) constraint_mx->makeDynamic());
          } else {
            throw _String("Dimensions of constraint matrix do not match current graph");
          }
        } // set node order
        else if (set_this_attribute == kBGMNodeOrder) {
          _Matrix     * order_mx  = (_Matrix *) _ProcessAnArgumentByType(*GetIthParameter(2UL), MATRIX, current_program, &dynamic_variable_manager);

          if (order_mx->check_dimension(1,num_nodes)) {

            _SimpleList order_list;
            order_mx->ConvertToSimpleList(order_list);

            bgm->SetNodeOrder ( &order_list );
          } else {
            throw _String("Order must be a row vector whose dimension matches the number of nodes in graph");
          }
        } else {
          throw (GetIthParameter (2UL)->Enquote() & " is not a valid parameter for BGM objects");
        }
      } // end BGM
      break;

      case HY_BL_SCFG:
      case HY_BL_LIKELIHOOD_FUNCTION: {


     if (object_type == HY_BL_SCFG && set_this_attribute == hy_env::kSCFGCorpus) {
          HBLObjectRef corpus_source = _ProcessAnArgumentByType (*GetIthParameter(2UL), MATRIX|STRING, current_program, &dynamic_variable_manager);
          if (corpus_source->ObjectClass () == STRING) {
            _List   single_string ( new _String (((_FString*)corpus_source)->get_str()));
            _Matrix wrapper (single_string);
            ((Scfg*)source_object)->SetStringCorpus (&wrapper);
          } else {
            _Matrix * matrix_corpus = (_Matrix*)corpus_source;
            if (matrix_corpus->IsAStringMatrix()) {
              ((Scfg*)source_object)->SetStringCorpus (matrix_corpus);
            } else {
              throw (set_this_attribute.Enquote() & " did not evaluate to a matrix of strings");
            }
          }
        } else {
          _LikelihoodFunction * lkf = (_LikelihoodFunction *) source_object;
          long parameter_index = _ProcessNumericArgumentWithExceptions (set_this_attribute ,current_program.nameSpacePrefix);
          if (lkf->GetIndependentVars().Map (parameter_index) < 0L) {
            throw (GetIthParameter(1)->Enquote() & " (=" & parameter_index & ") is not a valid parameter index");
           }
          lkf->SetIthIndependent (parameter_index,_ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix));
        }
      } // end HY_BL_SCFG; HY_BL_LIKELIHOOD_FUNCTION
        break;
      case HY_BL_DATASET:
      case HY_BL_DATASET_FILTER: {

        if (object_type == HY_BL_DATASET_FILTER) {
          ReleaseDataFilterLock (object_index);
        }

        long sequence_index  = _ProcessNumericArgumentWithExceptions (*GetIthParameter(1),current_program.nameSpacePrefix);
        _DataSet * ds = nil;
        if (object_type == HY_BL_DATASET) {
          ds = (_DataSet*) source_object;
        }
        else {
          _DataSetFilter *dsf = (_DataSetFilter*)source_object;
          ds = dsf->GetData ();
          sequence_index = dsf->theNodeMap.Map (sequence_index);
        }


        _String * sequence_name = new _String (_ProcessALiteralArgument (*GetIthParameter(2UL), current_program));

        if (! ds->SetSequenceName (sequence_index, sequence_name)) {
          delete sequence_name;
          throw  (GetIthParameter (1UL)->Enquote() & " (=" & sequence_index & ") is not a valid sequence index");
        }

      } // end HY_BL_DATASET; HY_BL_DATASET_FILTER
        break;

    } // end switch


  } catch (const _String& error) {
    return  _DefaultExceptionHandler (nil, error, current_program);
  }
  return true;
 }


  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleFprintf (_ExecutionList& current_program) {

  static const _String    kFprintfStdout               ("stdout"),
                          kFprintfMessagesLog          ("MESSAGE_LOG"),
                          kFprintfClearFile            ("CLEAR_FILE"),
                          kFprintfKeepOpen             ("KEEP_OPEN"),
                          kFprintfCloseFile            ("CLOSE_FILE"),
                          kFprintfSystemVariableDump   ("LIST_ALL_VARIABLES"),
                          kFprintfSelfDump             ("PRINT_SELF");


  static        _List       _open_file_handles_aux;
  static        _AVLListX   open_file_handles     (&_open_file_handles_aux);


  current_program.advance();

  bool     do_close                 = true,
  print_to_stdout          = false,
  skip_file_path_eval      = false,
  success                  = true;


  FILE*    destination_file = nil;

  try {

    _String  destination = *GetIthParameter(0UL);

    if (destination == kFprintfStdout) {
      _FString * redirect = (_FString*)hy_env::EnvVariableGet(hy_env::fprintf_redirect, STRING);
      if (redirect && redirect->has_data()) {
        destination         = redirect->get_str();
        if (destination == hy_env::kDevNull) {
          return true; // "print" to /dev/null
        } else {
          skip_file_path_eval = true;
        }
      }
      else {
        print_to_stdout = true;
      }
    } else {
        if (destination == kPromptForFilePlaceholder) {
            skip_file_path_eval = true;
        }
    }

    print_digit_specification = hy_env::EnvVariableGetDefaultNumber (hy_env::print_float_digits);

    if (!print_to_stdout) {
      if (destination == kFprintfMessagesLog) {
        if ((destination_file = hy_message_log_file) == nil) {
          return true; // requested print to MESSAGE_LOG, but it does not exist
                       // (e.g. remote MPI nodes, or running from a read only location
        }
      } else {
        if (skip_file_path_eval == false) {
          destination = _ProcessALiteralArgument (destination,current_program);
        }

        if (!ProcessFileName(destination, true,false,(hyPointer)current_program.nameSpacePrefix, false, &current_program)) {
          return false;
        }

        long open_handle  = open_file_handles.Find (&destination);

        do_close = open_handle < 0;

        if (!do_close) {
          destination_file = (FILE*)open_file_handles.GetXtra (open_handle);
        } else {
          if ((destination_file = doFileOpen (destination.get_str(), "a")) == nil)
            throw  (_String  ("Could not create/open output file at path ") & destination.Enquote() & ".");
        }
      }
    }

    for (unsigned long print_argument_idx = 1UL; print_argument_idx < parameter_count (); print_argument_idx++) {
      _String * current_argument = GetIthParameter(print_argument_idx);
      BaseRef   managed_object_to_print  = nil,
                dynamic_object_to_print = nil;

        // handle special cases first
      if (*current_argument == kFprintfClearFile) {
        if (!print_to_stdout && destination_file) {
          fclose (destination_file);
          destination_file = doFileOpen (destination.get_str(), "w");
            if (!do_close) {
              _String* destination_copy = new _String (destination);
            
              if (open_file_handles.UpdateValue(destination_copy, (long)destination_file, 1) >= 0) { // did not insert value
                DeleteObject (destination_copy);
              }
            }
          
        }
      } else if (*current_argument == kFprintfKeepOpen) {
        if (!print_to_stdout) {
          open_file_handles.Insert (new _String (destination), (long)destination_file, false, true);
          do_close = false;
        }
      } else if (*current_argument == kFprintfCloseFile) {
        open_file_handles.Delete (&destination, true);
        do_close = true;
      } else if (*current_argument == kFprintfSystemVariableDump ) {
        managed_object_to_print = &variableNames;
      } else if (*current_argument == kFprintfSelfDump) {
        managed_object_to_print = &current_program;
      } else {
          _String namespaced_id = AppendContainerName (*current_argument, current_program.nameSpacePrefix);

          // first check if the argument is a existing HBL object
          // TODO: this will go away in v3 when everything is in the same namespace


          if (!managed_object_to_print) { // not an existing variable
              long object_type = HY_BL_ANY, object_index;
              managed_object_to_print = //_GetHBLObjectByTypeMutable (namespaced_id, object_type, &object_index, false);
                _HYRetrieveBLObjectByNameMutable (namespaced_id, object_type,&object_index,false, false);
              if (managed_object_to_print) {
                  if (object_type == HY_BL_DATASET_FILTER) {
                    ReleaseDataFilterLock(object_index);
                  }
              } else {
                   dynamic_object_to_print = _ProcessAnArgumentByType (*current_argument, HY_ANY_OBJECT, current_program, nil);
              }
          }
      }

      BaseRef printables [2] = {managed_object_to_print, dynamic_object_to_print};
      for (BaseRef obj : printables) {
        if (obj) {
          if (!print_to_stdout) {
            obj->toFileStr (destination_file);
          } else {
             StringToConsole (_String((_String*)obj->toStr()));
          }
          break;
        }
      }
      DeleteObject(dynamic_object_to_print);

    } // end for

  }
  catch (const _String& error) {
    if (hy_env::EnvVariableTrue(hy_env::soft_fileio_exceptions)) {
         hy_env::EnvVariableSet(hy_env::last_fileio_exception, new _FString (error,false), false);
         success = true;
    } else {
        success =  _DefaultExceptionHandler (nil, error, current_program);
    }
  }

  if (destination_file && destination_file != hy_message_log_file && do_close) {
    fclose (destination_file);
  }

  return success;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleExecuteCommandsCases(_ExecutionList& current_program, bool do_load_from_file, bool do_load_library) {
    current_program.advance ();
    _String * source_code = nil;
    bool    pop_path = false;
    
    auto cleanup = [&] () -> void {
        if (pop_path) {
            PopFilePath();
        }
        DeleteObject(source_code);
    };
    
    _List dynamic_reference_manager;
    
    try {
        bool has_redirected_input = false,
        has_user_kwargs      = false;
        
        _List               _aux_argument_list;
        _AVLListXL          argument_list (&_aux_argument_list);
        _AssociativeList    * user_kwargs = nil;
        
        
        if (do_load_from_file) {
            _String file_path (*GetIthParameter(0UL));
            
            if (file_path == kPromptForFilePlaceholder ){
                ProcessFileName (file_path, false, false, (hyPointer)current_program.nameSpacePrefix);
            } else {
                file_path = _ProcessALiteralArgument(file_path, current_program);
            }
            _String original_path (file_path);

            FILE * source_file = nil;
            
            bool        reload          = hy_env::EnvVariableTrue(hy_env::always_reload_libraries);
            
            
            if (do_load_library) {
                bool has_extension    = file_path.FindBackwards (".",0,-1) != kNotFound;
                
                
                for (unsigned long p = 0; !source_file && p < _hy_standard_library_paths.countitems(); p++) {
                    for (unsigned long e = 0; !source_file && e < _hy_standard_library_extensions.countitems(); e++) {
                        _String try_path = *((_String*)_hy_standard_library_paths(p)) & file_path & *((_String*)_hy_standard_library_extensions(e));

                        ProcessFileName (try_path, false, false, (hyPointer)current_program.nameSpacePrefix, false, nil, false, true);
                        
                        if (loadedLibraryPaths.Find(&try_path) >= 0 && parameter_count() == 2UL && !reload) {
                            ReportWarning (_String("Already loaded ") & original_path.Enquote() & " from " & try_path);
                            return true;
                        }
                        if ((source_file = doFileOpen (try_path.get_str (), "rb"))) {
                            file_path = try_path;
                            break;
                        }
                        if (has_extension) {
                            break;
                        }
                    }
                }
            }
            
            if (source_file == nil) {
                ProcessFileName (file_path, false,false,(hyPointer)current_program.nameSpacePrefix);

                if (do_load_library && loadedLibraryPaths.Find(&file_path) >= 0 && parameter_count() == 2UL && !reload) {
                    ReportWarning (_String("Already loaded ") & original_path.Enquote() & " from " & file_path);
                    return true;
                }
                
                if ((source_file = doFileOpen (file_path.get_str (), "rb")) == nil) {
                    throw (_String("Could not read command file from '") &
                           original_path & "' (expanded to '" & file_path & "')");
                }
            }
            
            if (do_load_from_file && do_load_library && source_file) {
                ReportWarning (_String("Loaded ") & original_path.Enquote() & " from " & file_path.Enquote());
                loadedLibraryPaths.Insert (new _String (file_path),0,false,true);
            }
            
            source_code = new _String (source_file);
            
            if (fclose       (source_file) ) { // failed to fclose
                DeleteObject (source_code);
                throw (_String("Internal error: failed in a call to fclose ") & file_path.Enquote());
            }
            pop_path = true;
            PushFilePath (file_path);
        } else { // commands are not loaded from a file
            source_code = new _String (_ProcessALiteralArgument(*GetIthParameter(0UL), current_program));
        }
        
        
        if (!source_code || source_code->empty()) {
            throw _String("Empty/missing source code string");
        }
        
        if (!do_load_from_file && GetIthParameter(1UL)->nonempty()) {
            pop_path = true;
            PushFilePath (*GetIthParameter(1UL), false, false);
        }
        
        _String * use_this_namespace = nil;
        
        if (parameter_count() >= 3UL) { // stdin redirect (and/or name space prefix)
            _AssociativeList * input_arguments = nil;
            try {
                input_arguments =  (_AssociativeList *)_ProcessAnArgumentByType(*GetIthParameter(2UL), ASSOCIATIVE_LIST, current_program, &dynamic_reference_manager);
            } catch (const _String& err) {
                if (parameter_count() == 3UL) {
                    throw (err);
                }
            }
            
            
            if (input_arguments) {
                
                
                
                _List        *keys = input_arguments->GetKeys();
                dynamic_reference_manager < keys;
                keys->ForEach ([&] (BaseRef item, unsigned long) -> void {
                    _String * key = (_String*) item;
                    if (key) {
                        HBLObjectRef payload = input_arguments->GetByKey (*key, STRING | ASSOCIATIVE_LIST);
                        if (!payload) {
                            throw ((_String("All entries in the associative array used as input redirect argument to ExecuteCommands/ExecuteAFile must be strings or associative lists (for keyword arguments). The following key was not: ") & key->Enquote()));
                        }
                        if (key->BeginsWith ("--") && key->length() > 2) {
                            if (!user_kwargs) {
                                dynamic_reference_manager < (user_kwargs = new _AssociativeList);
                            }
                            user_kwargs->MStore(key->Cut (2, kStringEnd), payload, true);
                            has_user_kwargs = true;
                        } else {
                            argument_list.Insert (new _String (*key), (long)new _String (((_FString*)payload)->get_str()), false);
                            if (payload->ObjectClass() != STRING) {
                                throw ((_String("All entries in the associative array used as input redirect argument to ExecuteCommands/ExecuteAFile must be strings. The following key was not: ") & key->Enquote()));
                            }
                            has_redirected_input = true;
                        }
                    }
                });
                
                
                
                if (parameter_count() > 3UL) {
                    _String const namespace_for_code = _ProcessALiteralArgument (*GetIthParameter(3UL),current_program);
                    if (namespace_for_code.nonempty()) {
                        if (!namespace_for_code.IsValidIdentifier(fIDAllowCompound)) {
                            throw (_String("Invalid namespace ID in call to ExecuteCommands/ExecuteAFile: ") & GetIthParameter(3UL)->Enquote());
                        }
                        dynamic_reference_manager < (use_this_namespace = new _String (namespace_for_code));
                    }
                }
            }
        }
        
        if (parameter_count () < 4UL && current_program.nameSpacePrefix) {
            dynamic_reference_manager < (use_this_namespace = new _String (*current_program.nameSpacePrefix->GetName()));
        }
        
        if (source_code->BeginsWith ("#NEXUS")) {
            ReadDataSetFile (nil,1,source_code,nil,use_this_namespace);
        } else {
            bool result = false;
            
            _ExecutionList code (*source_code, use_this_namespace, false, &result);
            
            if (!result) {
                throw (_String("Encountered an error while parsing HBL"));
            } else {
                
                
                _AVLListXL * stash1 = nil;
                _List      * stash2 = nil,
                * stash_kw_tags = nil;
                
                _AssociativeList * stash_kw = nil;
                
                bool update_kw = false;
                
                if (has_redirected_input) {
                    code.stdinRedirectAux = &_aux_argument_list;
                    code.stdinRedirect = &argument_list;
                } else {
                    if (current_program.has_stdin_redirect()) {
                        stash1 = current_program.stdinRedirect;
                        stash2 = current_program.stdinRedirectAux;
                        current_program.stdinRedirect->AddAReference();
                        current_program.stdinRedirectAux->AddAReference();
                    }
                    code.stdinRedirect = current_program.stdinRedirect;
                    code.stdinRedirectAux = current_program.stdinRedirectAux;
                }
                
                
                bool ignore_ces_args = false;
                
                if (has_user_kwargs) {
                    code.SetKWArgs(user_kwargs);
                    ignore_ces_args = true;
                } else {
                    if (current_program.has_keyword_arguments()) {
                        code.kwarg_tags = stash_kw_tags = current_program.kwarg_tags;
                        code.kwargs = stash_kw = current_program.kwargs;
                        if (stash_kw_tags) current_program.kwarg_tags->AddAReference();
                        if (stash_kw) current_program.kwargs->AddAReference();
                        code.currentKwarg = current_program.currentKwarg;
                        update_kw = true;
                    }
                }
                
                
                
                if (!simpleParameters.empty() && code.TryToMakeSimple(true)) {
                    ReportWarning (_String ("Successfully compiled an execution list (possibly partially).\n") & _String ((_String*)code.toStr()) );
                    code.ExecuteSimple ();
                } else {
                    code.Execute(nil, ignore_ces_args);
                }
                
                if (stash1) {
                    stash1->RemoveAReference();
                    stash2->RemoveAReference();
                }
                
                if (stash_kw_tags) stash_kw_tags->RemoveAReference();
                if (stash_kw) stash_kw->RemoveAReference();
                
                code.stdinRedirectAux = nil;
                code.stdinRedirect    = nil;
                if (update_kw) {
                    code.kwarg_tags       = nil;
                    code.kwargs           = nil;
                    current_program.currentKwarg = code.currentKwarg;
                }
                
                if (code.result) {
                    DeleteObject (current_program.result);
                    current_program.result = code.result;
                    code.result = nil;
                }
            }
        }
    } catch (const _String& error) {
        cleanup ();
        return  _DefaultExceptionHandler (nil, error, current_program);
    }
    
    cleanup ();
    return true;
    
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleDoSQL (_ExecutionList& current_program) {

  static _SimpleList  sql_databases;

  static const _String kSQLOpen                 ("SQL_OPEN"),
                       kSQLClose                ("SQL_CLOSE");

  static auto sql_handler = [] (void* exL,int cc, char** rd, char** cn) -> int
    {
      static const _String  kSQLRowData              ("SQL_ROW_DATA"),
                            kSQLColumnNames          ("SQL_COLUMN_NAMES");

      _ExecutionList * caller = (_ExecutionList *)exL;

      if (!terminate_execution)
        if (caller && cc && !caller->empty()) {
          _List     row_data,
                    column_names;

          for (long cnt = 0; cnt < cc; cnt++) {
            if (rd[cnt]) {
              row_data.AppendNewInstance (new _String (rd[cnt]));
            } else {
              row_data.AppendNewInstance (new _String);
            }

            if (cn[cnt]) {
              column_names.AppendNewInstance (new _String (cn[cnt]));
            } else {
              column_names.AppendNewInstance (new _String);
            }
          }

          CheckReceptacleCommandIDException (&kSQLRowData,HY_HBL_COMMAND_DO_SQL,false)->SetValue (new _Matrix (row_data), false, true, NULL);
          CheckReceptacleCommandIDException (&kSQLColumnNames, HY_HBL_COMMAND_DO_SQL,false)->SetValue (new _Matrix (column_names), false,true, NULL);

          caller->Execute();

        }
       return 0;
    };


  _Variable * receptacle = nil;

  current_program.advance ();

  try {
#ifdef __HYPHY_NO_SQLITE__
    throw _String("SQLite commands can not be used in a HyPhy version built with the __HYPHY_NO_SQLITE__ flag");
#else

    if (*GetIthParameter(0UL) == kSQLOpen) { // open a DB from file path
      receptacle = _ValidateStorageVariable (current_program, 2UL);
      _String file_name = _ProcessALiteralArgument (*GetIthParameter(1UL), current_program);
      if (!ProcessFileName(file_name, true, false, (hyPointer)current_program.nameSpacePrefix,false,&current_program)) {
        return false;
      }

      sqlite3 *db = nil;
      int err_code = sqlite3_open (file_name.get_str(),&db);
      if (err_code == SQLITE_OK) {
        err_code = sqlite3_exec(db, "SELECT COUNT(*) FROM sqlite_master", sql_handler, nil, nil);
      }

      if (err_code != SQLITE_OK) {
        _String err_msg = sqlite3_errmsg(db);
        sqlite3_close(db);
        throw (err_msg);
      } else {
        long empty_slot = sql_databases.Find (0);
        if (empty_slot == kNotFound) {
          empty_slot = sql_databases.countitems();
          sql_databases << (long)db;
        } else {
          sql_databases.list_data[empty_slot] = (long)db;
        }

        sqlite3_busy_timeout (db, 5000);

        receptacle->SetValue (new _Constant (empty_slot), false,true, NULL);
      }
    } else {
      bool closing_db = *GetIthParameter(0UL) == kSQLClose;
      long db_index   = _ProcessNumericArgumentWithExceptions(*GetIthParameter(closing_db ? 2UL : 0UL), current_program.nameSpacePrefix);

      if (sql_databases.Map (db_index, 0L) == 0L) {
        throw (GetIthParameter(closing_db ? 2UL : 0UL)-> Enquote() & " is an invalid database index");
      }

      if (closing_db) {
        sqlite3_close ((sqlite3*)sql_databases.get(db_index));
        sql_databases [db_index] = 0L;
      } else {
          _String callback_code = _ProcessALiteralArgument(*GetIthParameter(2UL), current_program);

          _ExecutionList callback_script (callback_code,current_program.nameSpacePrefix?(current_program.nameSpacePrefix->GetName()):nil);

          if (!terminate_execution) {
            _String const sql_code = _ProcessALiteralArgument(*GetIthParameter(1UL), current_program);

            if (sqlite3_exec((sqlite3*)sql_databases.get(db_index), sql_code.get_str(), sql_handler, (hyPointer)&callback_script, nil) != SQLITE_OK) {
              throw _String(sqlite3_errmsg((sqlite3*)sql_databases(db_index)));
            }
          }
      }
    }
#endif

  } catch (const _String& error) {
    return  _DefaultExceptionHandler (receptacle, error, current_program);
  }

  return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleKeywordArgument (_ExecutionList& current_program) {
    current_program.advance ();
#ifdef __HYPHYMPI__
    // ignore keyword options for MPI nodes
    if (hy_mpi_node_rank > 0L) {
        return true;
    }
#endif
    
    try {
        _String keyword     = _ProcessALiteralArgument(*GetIthParameter(0UL), current_program),
                description = _ProcessALiteralArgument(*GetIthParameter(1UL), current_program),
                *default_value = nil;
        
        
        
        _List   reference_manager;
        if (parameter_count () > 2) {
            try {
                HBLObjectRef default_expression = _ProcessAnArgumentByType (*GetIthParameter(2UL), STRING|HY_UNDEFINED|NUMBER, current_program, nil);
                if (default_expression) {
                    reference_manager < default_expression;
                    if (default_expression->ObjectClass () == STRING) {
                        default_value = new _String (((_FString*)default_expression)->get_str());
                        reference_manager < default_value;
                    } else {
                        if (default_expression->ObjectClass () == NUMBER) {
                            default_value = new _String (((_Constant*)default_expression)->Value());
                            reference_manager < default_value;
                        } else {
                            // null values are handled separately
                        }
                    }
                }
            } catch (_String const& e){
                    throw _String ("Not a valid type for the default expression value");
            }
        }
        
        if (!current_program.kwarg_tags) {
            current_program.kwarg_tags = new _List;
        }
        
        _List * new_tagged_keyword = new _List (new _String (keyword), new _String (description));
        if (default_value) {
            (*new_tagged_keyword) << default_value;
        } else {
            if (parameter_count () > 2) {
                (*new_tagged_keyword) < (_String*)nil;
            }
        }

        if (parameter_count () > 3) {
            (*new_tagged_keyword) < new _String(_ProcessALiteralArgument(*GetIthParameter(3UL), current_program));
        }

        
        (*current_program.kwarg_tags) < new_tagged_keyword;
        
        //printf ("%s\n", _String ((_String*)current_program.kwarg_tags->toStr()).get_str());
        
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (nil, error, current_program);
    }
    
    return true;

}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetString (_ExecutionList& current_program) {

    auto make_fstring_pointer = [] (_String * s) -> _FString * {return new _FString (s);};
    auto make_fstring = [] (_String const s) -> _FString * {return new _FString (new _String (s));};

    static const _String kVersionString                   ("HYPHY_VERSION"),
                         kTimeStamp                       ("TIME_STAMP"),
                         kListLoadedLibraries             ("LIST_OF_LOADED_LIBRARIES");


    _Variable * receptacle = nil;

    current_program.advance ();

    try {
      receptacle = _ValidateStorageVariable (current_program);

      long         index1 = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2),current_program.nameSpacePrefix),
                   index2 = parameter_count() > 3 ? _ProcessNumericArgumentWithExceptions (*GetIthParameter(3),current_program.nameSpacePrefix) : -1L;


      // first, handle special, hardcoded cases

      HBLObjectRef return_value = nil;

      if (*GetIthParameter(1UL) == kVersionString) {
        if (index1 > 1.5) { // console header form
          return_value = make_fstring (_String ("HyPhy version ") & kHyPhyVersion);
        } else {
          if (index1 > 0.5) { // long form
            return_value = make_fstring(GetVersionString());
          } else {
            return_value = make_fstring (kHyPhyVersion);
          }
        }
      } else if (*GetIthParameter(1UL) == kTimeStamp) {
        return_value = make_fstring (GetTimeStamp (index1 < 0.5));
      }  else if (*GetIthParameter(1UL) == kListLoadedLibraries) {
        return_value = new _Matrix (loadedLibraryPaths.Keys());
      }

      if (!return_value) {

        //  next, handle cases like GetString (res, LikelihoodFunction, index)
          
        long       type_index = _HY_GetStringGlobalTypes.Find (GetIthParameter(1UL));

        if (type_index != kNotFound) {
          type_index = _HY_GetStringGlobalTypes.GetXtra (type_index);


          if (type_index != HY_BL_TREE) {
            _String * object_name = (_String*)GetObjectNameByType (type_index, index1);
            if (!object_name) {
              throw (_String ("There is no ") & GetIthParameter(1UL)->Enquote() & " object with index " & index1);
            }
            if (type_index == HY_BL_HBL_FUNCTION) {
              return_value =  &(
                                (*new _AssociativeList)
                              < (_associative_list_key_value){"ID", new _FString (*object_name) }
                              < (_associative_list_key_value){"Arguments", new _Matrix(GetBFFunctionArgumentList(index1)) }
                                );
            } else {
              return_value = make_fstring (*object_name);
            }

          } else {
            _String * tree_name = FetchMathObjectNameOfTypeByIndex (TREE, index1);
            if (!tree_name) {
              throw (_String ("There is no ") & GetIthParameter(1UL)->Enquote() & " object with index " & index1);
            }
            return_value = make_fstring(*tree_name);
          }

          receptacle->SetValue (return_value, false,true, NULL);
          return true;
        }

        // next, handle lookup of HBL objects
          
        const _String source_name   = AppendContainerName (*GetIthParameter(1), current_program.nameSpacePrefix);
        long          object_type = HY_BL_ANY,
                      object_index;
          
        BaseRefConst       source_object = nil;
        try {
          source_object = _GetHBLObjectByTypeMutable (source_name, object_type, &object_index);
        }
        catch (_String const & err) { // lookup failed
            return_value = nil;
            object_type = HY_BL_NOT_DEFINED;
        }

      switch (object_type) {
        case HY_BL_DATASET: {
          _DataSet const* data_set_object = (_DataSet const*)source_object;
          if (index1 >=0) { // get a specific sequence name
            if (index1 < data_set_object->NoOfSpecies()) {
              return_value = make_fstring (*(_String*)(data_set_object->GetNames().GetItem(index1)));
            } else {
              throw (_String (index1) & " exceeds the maximum index for the underlying DataSet object");
            }
          } else {
              return_value = new _Matrix (data_set_object->GetNames(), false);
          }
          break;
        }
        case HY_BL_DATASET_FILTER: {
          _DataSetFilter const* data_filter = (_DataSetFilter const*)source_object;
          ReleaseDataFilterLock (object_index);
          if (index1 >=0 ) { // get a specific sequence name
            if (index1 < data_filter->NumberSpecies()) {
              return_value = make_fstring (*(_String*)data_filter->GetData()->GetNames().GetItem(data_filter->theNodeMap.Element(index1)));
            } else {
              throw (_String (index1) & " exceeds the maximum index for the underlying DataSetFilter object");
            }
          } else {
            _List filter_seq_names;
            _List const * original_names = &data_filter->GetData()->GetNames();
            data_filter->theNodeMap.Each ([&] (long value, unsigned long) -> void {
              filter_seq_names << original_names->GetItem (value);
            });
            return_value = new _Matrix (filter_seq_names);

          }
          break;
        }

        case HY_BL_HBL_FUNCTION: {
          return_value  = &(
                            (*new _AssociativeList)
                            < (_associative_list_key_value){"ID", new _FString (*GetObjectNameByType (HY_BL_HBL_FUNCTION, object_index, false)) }
                            < (_associative_list_key_value){"Arguments", new _Matrix(GetBFFunctionArgumentList(object_index)) }
                            < (_associative_list_key_value){"Body", new _FString (GetBFFunctionBody(object_index).sourceText,false) }
                            );
          break;
        }

        case HY_BL_LIKELIHOOD_FUNCTION:
        case HY_BL_SCFG: {
          _LikelihoodFunction const *lf = (_LikelihoodFunction const*) source_object;
          if (index1 >= 0L) {
            if (index1<lf->GetIndependentVars().countitems()) {
              return_value = make_fstring (*LocateVar(lf->GetIndependentVars().GetElement(index1))->GetName());
            } else {
              if (index1 < lf->GetIndependentVars().countitems()+lf->GetDependentVars().countitems()) {
                return_value = make_fstring (*LocateVar(lf->GetDependentVars().GetElement(index1-lf->GetIndependentVars().countitems()))->GetName());
              } else {
                throw (_String (index1) & " exceeds the maximum index for the underlying LikelihoodFunction/SCFG object");
              }
            }
          } else {
            return_value = lf->CollectLFAttributes ();
            if (object_type == HY_BL_SCFG) {
              ((Scfg* const)lf)->AddSCFGInfo ((_AssociativeList*)return_value);
            }
          }
          break;
        }

        case HY_BL_MODEL: {
          if (index1 >= 0L) {
            if (index2 < 0L) { // get the name of index1 parameter
              long variable_index = PopulateAndSort ([&] (_AVLList & parameter_list) -> void  {
                  ScanModelForVariables (object_index, parameter_list, false, -1, false);
              }).Map (index1);
              if (variable_index >= 0) {
                return_value = make_fstring (*LocateVar(variable_index)->GetName());
              } else {
                throw (_String (index1) & " exceeds the maximum parameter index for the underlying Model object");
              }
            } else { // get the formula for cell (index1, index2)
              if (!IsModelOfExplicitForm (object_index)) {
                _Variable*      rate_matrix = (_Variable*)source_object;
                _Formula * cell = ((_Matrix*)rate_matrix->GetValue())->GetFormula (index1,index2);
                if (cell) {
                  return_value = make_fstring_pointer ((_String*)cell->toStr(kFormulaStringConversionNormal));
                } else {
                  throw (_String("Invalid rate matrix cell index"));
                }
              } else {
                throw (_String("Direct indexing of rate matrix cells is not supported for expression based (e.g. mixture) substitution models"));
              }
            }

          } else {
            _Variable   * rates, * freqs;
            bool         is_canonical;
            RetrieveModelComponents (object_index, rates, freqs, is_canonical);

            if (rates) {
              if (index1 == -1) { // branch length expression
                return_value = make_fstring_pointer (((_Matrix*)rates->GetValue())->BranchLengthExpression((_Matrix*)freqs->GetValue(),is_canonical));
              } else
              /*
               returns an AVL with keys
               "RATE_MATRIX" - the ID of the rate matrix
               "EQ_FREQS"    - the ID of eq. freq. vector
               "MULT_BY_FREQ" - a 0/1 flag to determine which format the matrix is in.
               */
              {
                return_value  = &(
                                  (*new _AssociativeList)
                                  < (_associative_list_key_value){"RATE_MATRIX",new _FString(*rates->GetName())}
                                  < (_associative_list_key_value){"EQ_FREQS",new _FString(*freqs->GetName()) }
                                  < (_associative_list_key_value){"MULT_BY_FREQ",new _Constant (is_canonical) }
                                  );
              }
            } else {
              throw _String("Failed to retrieve model rate matrix");
            }
          }
          break;
        }
        case HY_BL_BGM: {
            //ReportWarning(_String("In HandleGetString() for case HY_BL_BGM"));
          _BayesianGraphicalModel * this_bgm      = (_BayesianGraphicalModel *) source_object;

          switch (index1) {
            case HY_HBL_GET_STRING_BGM_SCORE: {   // return associative list containing node score cache
              _AssociativeList        * export_alist  = new _AssociativeList;

              if (this_bgm -> ExportCache (export_alist)) {
                return_value = export_alist;
              } else {
                DeleteObject (export_alist);
                throw _String("Failed to export node score cache for BGM");
              }

              break;
            }
            case HY_HBL_GET_STRING_BGM_SERIALIZE: {   // return associative list with network structure and parameters
              _StringBuffer *serialized_bgm = new _StringBuffer (1024L);
              this_bgm -> SerializeBGM (*serialized_bgm);
              return_value = new _FString (serialized_bgm);
              break;
            }
            default: {
              throw _String ("Unrecognized index ") & index1 & " for a BGM object";
            }
          }
          break;
        }
      }
    }

        // finally, try to look up a variable
      if (!return_value) {
        _Variable* var = FetchVar(LocateVarByName (AppendContainerName(*GetIthParameter(1UL), current_program.nameSpacePrefix)));

        if (var) {
            if (var->IsIndependent() && index1 != -3) {
              if (!var->has_been_set()) {
                  return_value = new _MathObject;
              } else {
                  return_value = make_fstring_pointer ((_String*) var->toStr());
              }
            } else {
              if (index1 == -1){
                _SimpleList variable_list = PopulateAndSort ([&] (_AVLList & parameter_list) -> void  {
                  var->ScanForVariables (parameter_list, true);
                });
                _AssociativeList   * var_list_by_kind = new _AssociativeList;

                _List split_vars;
                SplitVariableIDsIntoLocalAndGlobal (variable_list, split_vars);
                InsertVarIDsInList (var_list_by_kind, "Global", *(_SimpleList*)split_vars(0));
                InsertVarIDsInList (var_list_by_kind, "Local",  *(_SimpleList*)split_vars(1));
                return_value = var_list_by_kind;
              }
              else {  // formula string

                if (index1 == -3 || index1 == -4) {
                  _StringBuffer local, global;
                    
                  _SimpleList var_index;
                  if (index1 == -3 || var->IsIndependent())  {
                      var_index << var->get_index ();
                      if (var->IsIndependent()) {
                          //printf ("ExportIndVariables\n");
                        ExportIndVariables (global, local, &var_index);
                      } else {
                          //printf ("ExportDepVariables\n");
                        ExportDepVariables(global, local, &var_index);
                      }
                  } else {
                      _AVLList vl (&var_index);
                      var->ScanForVariables(vl, true);
                      _SimpleList ind_vars = var_index.Filter([] (long index, unsigned long) -> bool {return LocateVar(index)->IsIndependent();}),
                                  dep_vars = var_index.Filter([] (long index, unsigned long) -> bool {return !LocateVar(index)->IsIndependent();});
                      ExportIndVariables(global, local, &ind_vars);
                      ExportDepVariables(global, local, &dep_vars);
                  }
                  return_value = make_fstring_pointer ( & ((*new _StringBuffer (128L)) << global << local << '\n'));

                } else {
                  _Matrix * formula_matrix = (index2 >= 0 && var->ObjectClass() == MATRIX) ? (_Matrix*)var->GetValue () : nil;
                  if (formula_matrix) {
                    _Formula* cell = formula_matrix->GetFormula(index1, index2);
                    if (cell) {
                      return_value = make_fstring_pointer ((_String*) cell->toStr(kFormulaStringConversionNormal));
                    }
                  } else {
                    return_value = make_fstring_pointer ((_String*)var->GetFormulaString (kFormulaStringConversionNormal));
                  }
                }
              }
            }
          }
      }

      if (!return_value) {
        throw (_String("No viable object to obtain information from"));
      }

      receptacle->SetValue (return_value, false,true, NULL);


   }

    catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    
    return true;

}


  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleFscanf (_ExecutionList& current_program, bool is_sscanf) {
  
  const     static _String kFscanfStdin ("stdin");
  static    long   last_call_stream_position = 0L;
  
  current_program.advance();
  
  _List dynamic_reference_manager;
  
  try {
    
    long     current_stream_position = 0L,
             started_here_position   = 0L;
    
    unsigned long      has_rewind  = simpleParameters.Element(-1L) < 0 ? 1UL : 0UL;
    
    _String const * input_data = nil;
    
    _String source_name = *GetIthParameter(0UL);
    if (source_name == kFscanfStdin) {
      bool need_to_ask_user = true;
      if (current_program.has_stdin_redirect () || current_program.has_keyword_arguments()) {
          try {
            _FString * redirect = (_FString*)hy_env::EnvVariableGet(hy_env::fprintf_redirect, STRING);
            _String  * redirected = current_program.FetchFromStdinRedirect(nil, false, !(redirect && redirect->has_data()));
            input_data = redirected;
            // echo the input if there is no fprintf redirect in effect

            
            dynamic_reference_manager < redirected;
            need_to_ask_user = false;
          } catch (const _String& e) {
              if (e != kNoKWMatch) {
                  throw (e);
              }
          }
      }
        
      if (need_to_ask_user){ // read from stdin
        if (hy_env::EnvVariableTrue(hy_env::end_of_file) == false && source_name == hy_scanf_last_file_path)  {
          throw _String("Ran out of standard input");
        }
        _String * console_data = new _String (StringFromConsole());
        dynamic_reference_manager < console_data;
        input_data = console_data;
      }
    } else { // not stdin
      
      if (is_sscanf) {
        _FString * source_string = (_FString*)_ProcessAnArgumentByType(source_name,STRING,current_program,&dynamic_reference_manager);
        input_data = & source_string->get_str();
        if (hy_env::EnvVariableTrue(hy_env::end_of_file)) { // reset path
          hy_scanf_last_file_path = kEmptyString;
        }
        
        if (source_name != hy_scanf_last_file_path || has_rewind) { // new stream, or rewind on the current stream
          hy_scanf_last_file_path = source_name;
          current_stream_position = 0L;
          last_call_stream_position = 0L;
        } else {
          current_stream_position    = last_call_stream_position;
          started_here_position      = last_call_stream_position;
          if (last_call_stream_position >= input_data->length()) {
            // run out of input chars, set EOF and bail
            hy_env::EnvVariableSet(hy_env::end_of_file, new HY_CONSTANT_TRUE, false);
            return true;
          }
        }
      } else {
         _String file_path;
         if (source_name == kPromptForFilePlaceholder) {
             file_path = source_name;
         } else {
             file_path = ((_FString*)_ProcessAnArgumentByType(source_name,STRING,current_program,&dynamic_reference_manager))->get_str();
         }
         if (!ProcessFileName(file_path, true,false,(hyPointer)current_program.nameSpacePrefix, false, &current_program)) {
              return false;
         }
        
        FILE * input_stream = doFileOpen (file_path.get_str(), "rb");
        if (!input_stream) {
          throw     (file_path.Enquote() & " could not be opened for reading by fscanf. Path stack:\n\t" & GetPathStack("\n\t"));
        }
        if (hy_env::EnvVariableTrue(hy_env::end_of_file)) { // reset path
          hy_scanf_last_file_path = kEmptyString;
        }
        if (source_name != hy_scanf_last_file_path || has_rewind) { // new stream, or rewind on the current stream
          hy_scanf_last_file_path = source_name;
          last_call_stream_position = 0L;
        }

        fseek (input_stream,0,SEEK_END);
        current_stream_position    = ftell (input_stream);
        current_stream_position   -= last_call_stream_position;
        
        if (current_stream_position<=0) {
          hy_env::EnvVariableSet(hy_env::end_of_file, new HY_CONSTANT_TRUE, false);
          fclose(input_stream);
          return true;
        }
        
        rewind (input_stream);
        fseek  (input_stream, last_call_stream_position, SEEK_SET);
        _String * file_data = new _String (input_stream, current_stream_position);
        fclose (input_stream);
        dynamic_reference_manager < file_data;
        input_data = file_data;
        current_stream_position = 0L;
      }
    }
      
    unsigned long argument_index = 0UL,
                  upper_bound = has_rewind ? simpleParameters.countitems() - 1L : simpleParameters.countitems();
    
    while (argument_index < upper_bound && current_stream_position < input_data->length()) {
      _Variable * store_here = _ValidateStorageVariable (current_program, argument_index + 1UL);
     
      long   lookahead = current_stream_position;
      
      switch (simpleParameters.get(argument_index)) {
          
        case 0L: { // number
          _SimpleList numerical_match (input_data->RegExpMatch(hy_float_regex, lookahead));
          if (numerical_match.empty()) {
            throw (_String("Failed to read a number from the input stream") & input_data->Cut (lookahead, kStringEnd));
            break;
          }
          
          store_here->SetValue (new _Constant (input_data->Cut (numerical_match(0), numerical_match(1)).to_float ()), false,true, NULL);
          lookahead = input_data->FirstNonSpaceIndex(numerical_match (1) + 1, kStringEnd) ;
          }
          break;
          
        case 3L: { // string
          lookahead = 0L;
          bool  start_found=false;
          while (current_stream_position + lookahead < input_data->length ()) {
            char c = input_data->char_at (current_stream_position + lookahead);
            if (!start_found) {
              if (!isspace(c)) {
                current_stream_position += lookahead;
                start_found = true;
                lookahead = 0L;
              }
            } else if (c=='\n' || c=='\r' || c=='\t') {
              break;
            }
            lookahead++;
          }
          
          if (start_found) {
            store_here->SetValue (new _FString (new _String(*input_data,current_stream_position,current_stream_position+lookahead-1)),false,true, NULL);
          } else {
            store_here->SetValue (new _FString, false,true, NULL);
          }
          lookahead = current_stream_position + lookahead - 1L;
        }
        break;
        
        case 5L: { // Raw
          store_here->SetValue (new _FString (new _String (*input_data,current_stream_position,kStringEnd)), false,true, NULL);
          lookahead = input_data->length();
        }
        break;
        
        case 6L: { // Lines
          
          
          
          _String   line_block  (*input_data,current_stream_position,kStringEnd);
          
            // break the line block by any of the three platform line breaks
            // \r, \n or \r\n
          _List lines;
          long  last_break = 0L;
          
          auto add_buffer = [&] (long s, long e) -> void {
            if (e > s) {
              lines.AppendNewInstance(new _String (line_block, s, e-1));
            } else {
              lines.AppendNewInstance(new _String (kEmptyString));
            }
          };
          
          for (long i = 0UL; i < line_block.length(); i++) {
              char current_char = line_block.char_at(i);
              if (current_char == '\n') {
                add_buffer (last_break, i);
                last_break = i + 1L;
              } else {
                if (current_char == '\r') {
                  add_buffer (last_break, i);
                  if (line_block.char_at(i + 1L) == '\n') {
                    i ++;
                  }
                  last_break = i + 1L;
                }
              }
          }
          
          if (last_break < line_block.length()) {
            add_buffer (last_break, line_block.length ());
          }
        
          store_here->SetValue (new _Matrix (lines, false), false,true,NULL);
          lookahead = input_data->length();
          
        }
        
        break;
        
        default: {
            // TODO: 20170623 SLKP CHANGE: use expression extractor
          
          lookahead = input_data->ExtractEnclosedExpression(current_stream_position, (simpleParameters.list_data[argument_index]==2)?'(':'{', (simpleParameters.list_data[argument_index]==2)?')':'}', fExtractRespectQuote | fExtractRespectEscape);
          
          if (lookahead == kNotFound) {
            lookahead = input_data->length ();
            break;
          }
          
          _String object_data (*input_data, current_stream_position, lookahead);
          
          if (simpleParameters.list_data[argument_index] != 2) { // matrix
            _FormulaParsingContext def;
            store_here->SetValue (new _Matrix (object_data,simpleParameters.list_data[argument_index]==4, def), false,true,NULL);
          } else {
            _TheTree (*store_here->GetName(), object_data);
          }

        }
        break;
          
      } // end switch
      
      current_stream_position = lookahead + 1L;
      argument_index ++;
    }
    if (argument_index + has_rewind <simpleParameters.countitems()) {
      hy_env::EnvVariableSet(hy_env::end_of_file, new HY_CONSTANT_TRUE, false);
      throw (_String("Could not read all the parameters requested."));
    } else {
      hy_env::EnvVariableSet(hy_env::end_of_file, new HY_CONSTANT_FALSE, false);
    }
   
    last_call_stream_position += current_stream_position - started_here_position;

  } catch (const _String& error) {
      if (hy_env::EnvVariableTrue(hy_env::soft_fileio_exceptions)) {
          hy_env::EnvVariableSet(hy_env::last_fileio_exception, new _FString (error,false), false);
          return true;
      }
      return  _DefaultExceptionHandler (nil, error, current_program);
  }
  
  return true;
  
}

  //____________________________________________________________________________________

bool      _ElementaryCommand::HandleChoiceList (_ExecutionList& current_program) {
    
    auto   handle_exclusions = [] (long count, _SimpleList & excluded) -> const _SimpleList {
        _SimpleList lfids;
        lfids.Subtract(_SimpleList(count, 0L, 1L), excluded);
        return lfids;
    };
    
    static const _String kSkipNone ("SKIP_NONE"),
    kNoSkip ("NO_SKIP"),
    kLikelihoodFunctions ("LikelihoodFunction");
  
    static const unsigned long maximum_wrong_choices = 10UL;
  
    current_program.advance();
    
    _Variable * receptacle = nil;
    
    _List   local_dynamic_manager;
    
    
    try {
        
        receptacle = _ValidateStorageVariable (current_program);
        
        long    number_of_choices = _ProcessNumericArgumentWithExceptions (*GetIthParameter(2UL),current_program.nameSpacePrefix);
        _String dialog_title      = _ProcessALiteralArgument(*GetIthParameter(1UL), current_program),
                exclusions        = *GetIthParameter(3UL);
        
        _SimpleList selections,
        excluded;
        
        bool        variable_number = number_of_choices <= 0L,
                    do_markdown     = hy_env :: EnvVariableTrue(hy_env :: produce_markdown_output);
        
        
        if (exclusions != kSkipNone && exclusions != kNoSkip) {
            try {
                HBLObjectRef exlcusion_argument = _ProcessAnArgumentByType(exclusions, NUMBER | MATRIX, current_program, &local_dynamic_manager);
                if (exlcusion_argument->ObjectClass() == NUMBER) {
                    excluded << exlcusion_argument->Compute ()->Value();
                } else {
                    ((_Matrix*)exlcusion_argument)->ConvertToSimpleList (excluded);
                    excluded.Sort();
                }
            } catch (_String const & e) {
                //printf ("%s\n", e.get_str());
                // no exclusions, so do nothing
            }
        }
        
        _List  * available_choices = nil;
        
        if (simpleParameters.Element(0UL)) {
            // dynamically generated list of options
            _String const choices_parameter = *GetIthParameter(4UL);
            local_dynamic_manager < (available_choices = new _List);
            
            if (choices_parameter == kLikelihoodFunctions) {
                // the list consists of all defined likelihood function objects
                
                handle_exclusions (likeFuncNamesList.countitems(), excluded).Each([&] (long value, unsigned long) -> void {
                    if (likeFuncList.GetItem(value)) {
                        _String const * lf_name = (_String*) likeFuncNamesList (value);
                        (*available_choices) < new _List (new _String (*lf_name), new _String ( _String ("Likelihood Function ") & *lf_name & "."));
                    }
                });
            } else {
                // see if the argument is a reference to one of the standard HBL objects
                const _String source_name   = AppendContainerName (choices_parameter, current_program.nameSpacePrefix);
                
                long          object_type = HY_BL_DATASET_FILTER | HY_BL_DATASET | HY_BL_MODEL,
                object_index;
                
                
                try {
                    BaseRefConst       source_object = _GetHBLObjectByType (source_name, object_type, &object_index, &current_program);
                    // this wil also handle USE_LAST_MODEL
                    switch (object_type) {
                        case HY_BL_DATASET: {
                            _DataSet const *ds = (_DataSet const*) source_object;
                            handle_exclusions (ds->NoOfSpecies(), excluded).Each (
                                                                                  [&] (long value, unsigned long) -> void {
                                                                                      _String const * sequence_name = ds->GetSequenceName(value);
                                                                                      (*available_choices) < new _List (new _String (*sequence_name), new _String ( _String ("Taxon ") & (value + 1L) & sequence_name->Enquote('(', ')') & "."));
                                                                                  }
                                                                                  );
                            break;
                        }
                        case HY_BL_DATASET_FILTER: {
                            _DataSetFilter const *df = (_DataSetFilter const*) source_object;
                            handle_exclusions (df->NumberSpecies(), excluded).Each (
                                                                                    [&] (long value, unsigned long) -> void {
                                                                                        _String const * sequence_name = df->GetSequenceName(value);
                                                                                        (*available_choices) < new _List (new _String (*sequence_name), new _String ( _String ("Taxon ") & (value + 1L) & sequence_name->Enquote('(', ')') & "."));
                                                                                    }
                                                                                    );
                            break;
                        }
                        case HY_BL_MODEL: {
                            (*available_choices) < new _List (new _String("All Parameters"), new _String("All local model parameters are constrained"));
                            
                            _SimpleList model_indices = PopulateAndSort ([&] (_AVLList & parameter_list) -> void  {
                                if (IsModelOfExplicitForm (object_index)) {
                                    ((_Formula const*)source_object)->ScanFForVariables(parameter_list,false);
                                } else {
                                    ((_Variable const*)source_object)->ScanForVariables(parameter_list,false);
                                }
                            });
                            handle_exclusions (model_indices.countitems(), excluded).Each (
                                                                                           [&] (long value, unsigned long) -> void {
                                                                                               _String const * parameter_name = LocateVar(model_indices.get(value))->GetName();
                                                                                               (*available_choices) < new _List (new _String (*parameter_name),
                                                                                                                                 new _String (_String ("Constrain parameter ") & *parameter_name & '.'));
                                                                                           }
                                                                                           );
                        }
                            break;
                    }
                } catch (_String const & e) {
                    // not an object
                    try {
                        _Matrix * target_variable = (_Matrix*)_ProcessAnArgumentByType(choices_parameter, MATRIX, current_program, &local_dynamic_manager);
                        if (!target_variable->IsAStringMatrix()) {
                            throw (choices_parameter.Enquote() & " is not a matrix of strings");
                        }
                        if (target_variable->GetVDim() != 2) {
                            throw (choices_parameter.Enquote() & " is not a matrix with two columns");
                        }
                        
                        _List choices;
                        target_variable->FillInList(choices, false);
                        //choices.bumpNInst();
                        
                        handle_exclusions (target_variable->GetHDim(), excluded).Each (
                             [&] (long value, unsigned long ) -> void {
                                 _String const * parameter_name = LocateVar(value)->GetName();
                                 BaseRef key = choices.GetItem(value << 1),
                                         description = choices.GetItem(1L + (value << 1));
                                 key->AddAReference(); description->AddAReference();
                                 (*available_choices) < new _List (key,description);
                             }
                             );
                                                                                       
                        
                        /*for (unsigned long k = 0UL; k < choices.countitems(); k+=2) {
                            (*available_choices) < new _List (choices.GetItem(k), choices.GetItem(k+1));
                        }*/
                    }   catch (_String const& e2) {
                        throw (choices_parameter.Enquote() & " is not a supported object/literal for the list of choices");
                    }
                    
                }
                
            }
            
        } else {
            available_choices = (_List*)parameters.GetItem(4UL);
        }
        
        if (available_choices->empty()) {
            throw (_String("The list of choices is empty"));
        }
        
        if ((long) available_choices->countitems() < number_of_choices) {
            throw (_String ("The list of choices has " ) & (long) available_choices->countitems() & " elements, but " & number_of_choices & " choices are required");
        }
        
        auto validate_choice = [&] (_String const& user_choice) -> long {
            return available_choices->FindOnCondition([&] (BaseRefConst item, unsigned long) -> bool {
                return user_choice == *(_String*)((_List*)item)->GetItem (0);
            });
        };
        
        long required = variable_number ? available_choices->countitems() : number_of_choices;
        bool need_to_prompt_user = true;

        if (current_program.has_stdin_redirect() || current_program.has_keyword_arguments()) {
            
            auto report_key_error = [&] (const _String & bad_key) -> void {
                _String choice_list_echo = available_choices->Map ([] (BaseRef object, unsigned long index) -> BaseRef {
                    BaseRef title = ((_List*)object)->GetItem (0);
                    title->AddAReference();
                    return title;
                }, 0, MIN (32, available_choices->countitems())).Join(", ");
                if (available_choices->countitems() > 32) {
                    choice_list_echo = choice_list_echo & " and " & _String(available_choices->countitems()-32) & " additional options";
                }
                
                throw (bad_key.Enquote() & " is not a valid choice passed to '" & dialog_title & "' ChoiceList using redirected stdin input or keyword arguments. Valid choices are\n\t") & choice_list_echo & "\n";
            };
            
            
            while (selections.countitems() < required) {
                _String user_choice;
                try {
                    _FString * redirect = (_FString*)hy_env::EnvVariableGet(hy_env::fprintf_redirect, STRING);
                    user_choice = current_program.FetchFromStdinRedirect(&dialog_title, required > 1, !(redirect && redirect->has_data())); // handle multiple selections
                } catch (const _String& e) {
                    if (e == kNoKWMatch) {
                        break;
                    } else {
                        throw (e);
                    }
                } catch (HBLObjectRef multiple_choice) {
                    if (multiple_choice->ObjectClass() == ASSOCIATIVE_LIST) {
                        _List option_list;
                        ((_AssociativeList*)multiple_choice)->FillInList (option_list);
                        option_list.ForEach([&] (BaseRef item, unsigned long) -> void {
                            long    match_found = validate_choice (_String ((_String*)item->toStr()));
                            if (match_found == kNotFound)  {
                                report_key_error (_String ((_String*)item->toStr()));
                            } else {
                                selections << match_found;
                            }
                        });
                        DeleteObject(multiple_choice);
                    }
                }
                
                if (variable_number && user_choice.empty()) {
                    break;
                }
                
                long    match_found = validate_choice (user_choice);
                if (match_found == kNotFound) {
                    report_key_error (user_choice);
                }
                selections << match_found;
            }
            
            need_to_prompt_user = selections.countitems() != required && !variable_number || variable_number && selections.empty();
        }
        
        if (need_to_prompt_user) {
#ifdef  __HEADLESS__
            throw (_String("Unhandled request for data from standard input (headless HyPhy)"));
#endif
            if (do_markdown) {
                printf ("\n\n####%s\n", dialog_title.get_str());
            } else {
                const _String spacer (_String('-'), dialog_title.length());
                printf ("\n\n\t\t\t+%s+\n\t\t\t|%s|\n\t\t\t+%s+\n\n",
                         spacer.get_str(),
                         dialog_title.get_str(),
                         spacer.get_str());
            }
          
            unsigned long wrong_selections = 0UL;
          
            _AVLList      selection_tree   (&selections);
          
          
            while (selections.countitems() < required) {
              
              for (unsigned long option = 0UL; option < available_choices->countitems(); option ++) {
                  if (selection_tree.FindLong (option) == kNotFound) {
                    printf (do_markdown ? "\n%ld. [**%s**] %s" : "\n\t(%ld):[%s] %s", option+1UL, ((_String*)available_choices->GetItem (option, 0))->get_str(), ((_String*)available_choices->GetItem (option, 1))->get_str());
                  }
                
               }
              
               if (variable_number) {
                  printf ("\n\n%sPlease choose option %ld, enter d to complete selection, enter q to cancel:", do_markdown ? ">" : "",selections.countitems() + 1UL);
               } else {
                if (number_of_choices > 1L) {
                  printf ("\n\n%sPlease choose option %ld of %ld (or press q to cancel selection):",do_markdown ? ">" : "", selections.countitems() + 1UL ,number_of_choices);
                } else {
                  printf ("\n\n%sPlease choose an option (or press q to cancel selection):", do_markdown ? ">" : "");
                 }
               }
              
               _String user_choice (StringFromConsole());
               if (user_choice.length () == 1UL && toupper(user_choice.char_at(0L)) == 'Q') {
                 selections.Clear();
                 selections << -1L;
                 break;
               } else {
                 if (variable_number && user_choice.length () == 1UL && toupper(user_choice.char_at(0L)) == 'D') {
                   if (selections.empty()) {
                     selections << -1L;
                   }
                   break;
                 }
               }
              
               long integer_choice = user_choice.to_long();
              
               if   (integer_choice > 0L && integer_choice <= available_choices->countitems()) {
                 selection_tree.Insert((BaseRef)(integer_choice-1));
               } else {
                 ++ wrong_selections;
               }
             
               if (wrong_selections > maximum_wrong_choices) {
                throw _String ("Too many invalid imputs");
               }
            }
          
            selection_tree.ReorderList();
          
          
        }
      
        if (selections.empty () || selections.GetElement(0UL) < 0) {
            // failed selection
            hy_env::EnvVariableSet(hy_env::selection_strings, new HY_NULL_RETURN, false);
            if (number_of_choices == 1L) {
                receptacle->SetValue (new _Constant (-1.), false,true,NULL);
            } else {
                receptacle->SetValue (new _Matrix (_SimpleList (1UL,-1L,0)), false,true,NULL);
            }
            terminate_execution = true;
        } else {
            if (variable_number || number_of_choices > 1L) {
                _SimpleList corrected_for_exclusions;
                _List       chosen_strings;
                for (unsigned long i = 0UL; i < selections.countitems(); i++) {
                    corrected_for_exclusions << excluded.SkipCorrect(selections.Element(i));
                    chosen_strings < new _FString (*(_String*)available_choices->GetItem(selections.Element(i), 0),false);
                }
                receptacle->SetValue (new _Matrix (corrected_for_exclusions), false, true, NULL);
            } else {
                receptacle->SetValue (new _Constant (excluded.SkipCorrect (selections.Element(0UL))), false, true, NULL);
                hy_env::EnvVariableSet(hy_env::selection_strings, new _FString (*(_String*)available_choices->GetItem(selections.Element(0UL), 0),false), false);
            }
        }
        
        
    } catch (const _String& error) {
        return  _DefaultExceptionHandler (receptacle, error, current_program);
    }
    
    return true;
    
}
  //____________________________________________________________________________________
  // REQUIRES CODE REVIEW FROM THIS POINT ON
  //____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase38 (_ExecutionList& chain, bool sample) {
  chain.currentCommand++;
  
  _List local_object_manager;
  
  _String *likef          =  GetIthParameter(1);
  
  _String name2lookup = AppendContainerName(*likef,chain.nameSpacePrefix);
  long    objectID    = FindLikeFuncName (name2lookup);
  try {
    if (objectID >= 0) {
      _DataSet     * ds               = new _DataSet;
      _String      * dsName           = new _String (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix));
      _LikelihoodFunction *lf         = ((_LikelihoodFunction*)likeFuncList(objectID));
      
      local_object_manager < ds < dsName;
      
      _Matrix * partitionList         = nil;
      if (parameters.lLength>2) {
        _String  secondArg = *GetIthParameter(2);
        local_object_manager < (partitionList = (_Matrix*)ProcessAnArgumentByType (&secondArg, chain.nameSpacePrefix, MATRIX));
      }
      _SimpleList                     partsToDo;
      
      if (lf->ProcessPartitionList(partsToDo, partitionList)) {
        lf->ReconstructAncestors(*ds, partsToDo, *dsName,  sample, simpleParameters.Find(-1) >= 0, simpleParameters.Find(-2) >= 0 );
      }
      
      ds->AddAReference();
      StoreADataSet  (ds, dsName);
    
    } else {
      objectID    =   FindSCFGName       (name2lookup);
      if (objectID>=0)
      /* reconstruct best parse tree for corpus using SCFG */
      {
        CheckReceptacleAndStore (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix)," ReconstructAncestors (SCFG)", true, new _FString( ((Scfg*)scfgList (objectID))->BestParseTree() ), false);
      } else {
        throw  (_String("Likelihood Function/SCFG") & *likef & _String(" has not been initialized"));
      }
    }
  } catch (const _String& error) {
    _DefaultExceptionHandler (nil, error, chain);
  }
}


  //____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase31 (_ExecutionList& chain) {
    // 20100312 SLKP: added matrix-expression based model
    // definitions
  chain.advance();
  
  try {
      // first check to see if matrix parameters here are valid
    
    bool     usingLastDefMatrix = false,
    doExpressionBased  = false;
    
    _Formula *isExpressionBased  = nil;
    
    _String* parameterName,
    errMsg,
    arg0 = chain.AddNameSpaceToID(*(_String*)parameters(0));
    
    long     f,
    f2=-1L,
    matrixDim,
    f3,
    multFreqs = 1;
    
    
    
    if (parameters.lLength>3) {
      parameterName = (_String*)parameters.list_data[3];
      if (parameterName->Equal(explicitFormMExp)) {
        doExpressionBased = true;
        multFreqs         = 0;
      } else {
        multFreqs = ProcessNumericArgument (parameterName,chain.nameSpacePrefix);
      }
    }
    
    _Matrix*  checkMatrix = nil;
    
    parameterName = (_String*)parameters.list_data[1];
    
    if (parameterName->Equal (useLastDefinedMatrix)) {
      if (lastMatrixDeclared<0) {
        throw _String ("First Call to Model. USE_LAST_DEFINED_MATRIX is meaningless.");
      }
      f3 = lastMatrixDeclared;
      f  = modelMatrixIndices[f3];
      usingLastDefMatrix = true;
    }
    
    
    if (doExpressionBased) {
      _String matrixExpression (ProcessLiteralArgument((_String*)parameters.list_data[1],chain.nameSpacePrefix)),
      defErrMsg = _String ("The expression for the explicit matrix exponential passed to Model must be a valid matrix-valued HyPhy formula that is not an assignment") & ':' & matrixExpression;
        // try to parse the expression, confirm that it is a square  matrix,
        // and that it is a valid transition matrix
      isExpressionBased = new _Formula;
      _FormulaParsingContext fpc (nil, chain.nameSpacePrefix);
      _StringBuffer  trimmed_expression;
      _ElementaryCommand::FindNextCommand (matrixExpression,trimmed_expression);
      matrixExpression = trimmed_expression;
      long parseCode = Parse(isExpressionBased,matrixExpression,fpc, nil);
      if (parseCode != HY_FORMULA_EXPRESSION || isExpressionBased->ObjectClass()!= MATRIX ) {
        throw (defErrMsg & " parse code = " & parseCode & " " & (parseCode == HY_FORMULA_EXPRESSION ? (_String(", object type code ") & _String((long) isExpressionBased->ObjectClass())) : kEmptyString ));
      }
      
        //for (unsigned long k = 0; k < isExpressionBased
      
      checkMatrix = (_Matrix*)isExpressionBased->Compute();
      
      
    } else {
      parameterName = (_String*)parameters.list_data[1];
      
      _String augName (chain.AddNameSpaceToID(*parameterName));
      f = LocateVarByName (augName);
      
      if (f<0) {
        throw (*parameterName & " has not been defined prior to the call to Model = ...");
      }
      
      _Variable* checkVar = usingLastDefMatrix?LocateVar(f):FetchVar (f);
      if (checkVar->ObjectClass()!=MATRIX) {
        throw (*parameterName & " must refer to a matrix in the call to Model = ...");
      }
      checkMatrix = (_Matrix*)checkVar->GetValue();
    }
    
    
    
      // so far so good
    matrixDim = checkMatrix->GetHDim();
    if ( matrixDim!=checkMatrix->GetVDim() || matrixDim<2 ) {
      throw (*parameterName & " must be a square matrix of dimension>=2 in the call to Model = ...");
    }
    
    parameterName = (_String*)parameters.list_data[2]; // this is the frequency matrix (if there is one!)
    _String         freqNameTag (chain.AddNameSpaceToID(*parameterName));
    
    f2 = LocateVarByName (freqNameTag);
    if (f2<0) {
      throw(*parameterName & " has not been defined prior to the call to Model = ...");
    }
    _Variable * checkVar = FetchVar (f2);
    if (checkVar->ObjectClass()!=MATRIX) {
      throw (*parameterName & " must refer to a column/row vector in the call to Model = ...");
     }
    
    checkMatrix = (_Matrix*)checkVar->GetValue();
    
    if (checkMatrix->GetVDim()==1UL) {
      if (checkMatrix->GetHDim()!=matrixDim) {
        throw (*parameterName & " must be a column vector of the same dimension as the model matrix in the call to Model = ...");
      }
    } else if (checkMatrix->GetHDim()==1UL) {
      if (checkMatrix->GetVDim()!=matrixDim) {
        throw ( *parameterName & " must be a row vector of the same dimension as the model matrix in the call to Model = ...");
      }
      errMsg = *parameterName & " has been transposed to the default column vector setting ";
      checkMatrix->Transpose();
      ReportWarning (errMsg);
    } else {
      throw (*parameterName & " must refer to a column/row vector in the call to Model = ...");
    }
    
    if (usingLastDefMatrix) {
      if (modelFrequenciesIndices[f3]<0) {
        f2 = -f2-1;
      }
    } else if (multFreqs == 0) { // optional flag present
      f2 = -f2-1;
    }
    
    long existingIndex = modelNames.FindObject(&arg0);
    
    if (existingIndex == -1) { // name not found
      lastMatrixDeclared = modelNames.FindObject (&kEmptyString);
      
      if (lastMatrixDeclared>=0) {
        modelNames.Replace (lastMatrixDeclared,&arg0,true);
        modelTypeList.list_data[lastMatrixDeclared] = isExpressionBased?matrixDim:0;
        if (isExpressionBased) {
          modelMatrixIndices.list_data[lastMatrixDeclared] = (long)isExpressionBased;
        } else {
          modelMatrixIndices.list_data[lastMatrixDeclared] = (usingLastDefMatrix?f:variableNames.GetXtra(f));
        }
        
        if (f2>=0) {
          modelFrequenciesIndices.list_data[lastMatrixDeclared] = variableNames.GetXtra(f2);
        } else {
          modelFrequenciesIndices.list_data[lastMatrixDeclared] = -variableNames.GetXtra(-f2-1)-1;
        }
      } else {
        modelNames && & arg0;
        modelTypeList << (isExpressionBased?matrixDim:0);
        if (isExpressionBased) {
          modelMatrixIndices << (long)isExpressionBased;
        } else {
          modelMatrixIndices << (usingLastDefMatrix?f:variableNames.GetXtra(f));
        }
        if (f2>=0) {
          modelFrequenciesIndices << variableNames.GetXtra(f2);
        } else {
          modelFrequenciesIndices << -variableNames.GetXtra(-f2-1)-1;
        }
        lastMatrixDeclared = modelNames.lLength-1;
      }
    } else {
      modelNames.Replace(existingIndex,&arg0,true);
      if (modelTypeList.list_data[existingIndex]) {
        delete ((_Formula*)modelMatrixIndices[existingIndex]);
      }
      
      modelTypeList.list_data[existingIndex] = isExpressionBased?matrixDim:0;
      if (isExpressionBased) {
        modelMatrixIndices[existingIndex] = (long)isExpressionBased;
      } else {
        modelMatrixIndices[existingIndex] = usingLastDefMatrix?f:variableNames.GetXtra(f);
      }
      
      
      if (f2>=0) {
        modelFrequenciesIndices[existingIndex] = variableNames.GetXtra(f2);
      } else {
        modelFrequenciesIndices[existingIndex] = -variableNames.GetXtra(-f2-1)-1;
      }
      
      lastMatrixDeclared = existingIndex;
    }
  } catch (const _String& error) {
    _DefaultExceptionHandler (nil, error, chain);
  }
}





