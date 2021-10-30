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

#include "defines.h"
#include "global_object_lists.h"
#include "hbl_env.h"
#include "function_templates.h"
#include "variable.h"
#include "batchlan.h"
#include "likefunc.h"
#include "global_things.h"



/** Legacy declarations */

extern _List batchLanguageFunctionNames;


_List        _data_filter_aux;
_AVLListXL   _data_filters (&_data_filter_aux);

_SimpleList  _data_filter_locks_aux;
_AVLList     _data_filter_locks (&_data_filter_locks_aux);

_SimpleList  _data_filter_listeners_aux;
_AVLListXL   _data_filter_listeners (&_data_filter_listeners_aux);

namespace hyphy_global_objects {
  /**
   notification types enum
   */
  
  enum kNotificationType {
    kNotificationTypeChange,
    kNotificationTypeDelete
  };
  
  /**
   data filter objects -- internal
   */
  
    _List _batchLanugageFunctionNamesIndexed;
   _AVLListX batchLanguageFunctionNamesIndexed (&_batchLanugageFunctionNamesIndexed);

  
  void  _SetDataFilterParameters (_String const& name, _DataSetFilter const& filter) {
    setParameter (WrapInNamespace ("species", &name), filter.NumberSpecies());
    setParameter (WrapInNamespace ("sites", &name), filter.GetSiteCountInUnits());
    
    hyFloat size_cutoff = hy_env::EnvVariableGetNumber(hy_env::dataset_save_memory_size);
    if (filter.GetSiteCount() < size_cutoff) {
      setParameter(WrapInNamespace("site_map", &name), new _Matrix (filter.theOriginalOrder), nil, false);
      setParameter(WrapInNamespace("site_freqs", &name), new _Matrix (filter.theFrequencies), nil, false);
    }
    
    if (filter.NumberSpecies() < size_cutoff) {
      setParameter(WrapInNamespace("sequence_map", &name), new _Matrix (filter.theNodeMap), nil, false);
    }
  }
  
  void      _KillDataFilterParameters (_String const & filter_name) {
    _List kill_these_arguments;
    kill_these_arguments < "species" < "sites" < "site_map" < "site_freqs" < "sequence_map";
    
    for (unsigned long i = 0; i < kill_these_arguments.countitems(); i++) {
      DeleteVariable(WrapInNamespace(*(_String*)kill_these_arguments.GetItem(i), &filter_name));
    }
  }
  
  void    _NotifyDataFilterListeners (const long index, kNotificationType event_type) {
    _List * listeners = (_List*)_data_filter_listeners.GetDataByKey( (BaseRef) index);
    if (listeners) {
      
     _List buffered_updates;
        
      for (unsigned long k = 0UL; k < listeners->lLength; k++) {
        BaseRef this_listener = listeners->GetItem(k);
        
        if (_LikelihoodFunction* lf = dynamic_cast<_LikelihoodFunction*> (this_listener)) {
            //StringToConsole(_String("_NotifyDataFilterListeners ") & index & " " & (long)lf & "\n");
       
          if (event_type == kNotificationTypeChange) {
            buffered_updates << lf;
            //lf->Rebuild();
               /* 20170328 SLKP: this COULD MODIFY the 'listeners' object, hence the buffering */
          } else if (event_type == kNotificationTypeDelete) {
            hy_global::HandleApplicationError (_String("Attempted to delete a data set filter ") & GetFilterName(index)->Enquote() & " which is referenced by the likelihood function " &  GetObjectNameByType (HY_BL_LIKELIHOOD_FUNCTION, FindLikeFuncIndex (lf), false)->Enquote());
          }
        }
      }
        
        for (unsigned long k = 0UL; k < buffered_updates.lLength; k++) {
            ((_LikelihoodFunction*)buffered_updates.GetItem (k)) -> Rebuild();
        }
    }
  }
  
  
  /**
   generic hidden functions
   */
  
  bool   _IsObjectLocked (long index, long object_class) {
    
    _AVLList* source = nil;
    
    switch (object_class) {
      case HY_BL_DATASET_FILTER:
        source = &_data_filter_locks;
        break;
    }
    
    if (source) {
      return source->FindLong (index) >= 0L;
    }
    
    return false;
  }
  
  bool   _AcquireObjectLock (long index, long object_class) {
    
    _AVLList* source = nil;
    
    switch (object_class) {
      case HY_BL_DATASET_FILTER:
        source = &_data_filter_locks;
        break;
    }
    
    if (source) {
      long added_as = source->Insert((BaseRef)index);
      return added_as >= 0;
    }
    
    return false;
  }
  
  bool   _RemoveObjectLock (long index, long object_class) {
    
    _AVLList* source = nil;
    
    switch (object_class) {
      case HY_BL_DATASET_FILTER:
        source = &_data_filter_locks;
        break;
    }
    
    if (source) {
      source->Delete ((BaseRef)index);
      return true;
    }
    
    return false;
  }
  
  /**
   public facing functions for data filter objects
   */
  
  const   _DataSetFilter * GetDataFilter (long index) {
    if (_data_filters.IsValidIndex (index)) {
      return (const _DataSetFilter*)_data_filters.GetXtra (index);
    }
    return nil;
  }
  
  const   _DataSetFilter * GetDataFilter (_String const& name ) {
    return GetDataFilter(FindDataFilter (name));
  }
  
  _DataSetFilter * ExclusiveLockDataFilter (long index) {
    if (_data_filters.IsValidIndex (index)) {
      if (_AcquireObjectLock(index, HY_BL_DATASET_FILTER)) {
        return (_DataSetFilter*)_data_filters.GetXtra (index);
      }
    }
    return nil;
  }
  
  _DataSetFilter * ExclusiveLockDataFilter (_String const& name) {
    return ExclusiveLockDataFilter(FindDataFilter(name));
  }
  
  bool    ReleaseDataFilterLock (long index) {
    if (_data_filters.IsValidIndex (index)) {
      return _RemoveObjectLock(index, HY_BL_DATASET_FILTER);
    }
    return false;
  }
  
  bool    ReleaseDataFilterLock (_String const& name) {
    return ReleaseDataFilterLock(FindDataFilter(name));
  }
  
  
  bool    UnregisterChangeListenerForDataFilter (long const index, BaseRef listener) {
    if (_data_filters.IsValidIndex (index)) {
      
      if (dynamic_cast <_LikelihoodFunction*> (listener)) {
        _List * current_listeners = (_List *)_data_filter_listeners.GetDataByKey ((BaseRef)index);
        if (current_listeners) {
          long listener_index = current_listeners->_SimpleList::Find((long)listener);
          if (listener_index >= 0) {
            //StringToConsole(_String("UnregisterChangeListenerForDataFilter ") & index & " " & (long)listener & "\n");
            current_listeners->Delete (listener_index);
          }
         return true;
        }
      }
       hy_global::HandleApplicationError (_String("Not a supported listener type in call to ") & _String (__PRETTY_FUNCTION__));
    }
    
    return false;
  }
  
  
  bool    RegisterChangeListenerForDataFilter (long const index, BaseRef listener) {
      
    //StringToConsole(_String("RegisterChangeListenerForDataFilter ") & index & " " & (long)listener & "\n");
      
    if (_data_filters.IsValidIndex (index)) {
      
      if (dynamic_cast <_LikelihoodFunction*> (listener)) {
        _List * current_listeners = (_List *)_data_filter_listeners.GetDataByKey ((BaseRef)index);
        if (!current_listeners) {
          current_listeners = new _List;
          _data_filter_listeners.Insert ((BaseRef)index, (long)current_listeners, false, false);
        }
        
        
        if (current_listeners->_SimpleList::Find((long)listener) < 0L) {
          (*current_listeners) << listener;
        }
        return true;
        
      }
      
      
       hy_global::HandleApplicationError (_String("Not a supported listener type in call to ") & _String (__PRETTY_FUNCTION__));
      
      
    }
    
    return false;
  }
  
  
  bool    RegisterChangeListenerForDataFilter (_String const& name, BaseRef listener) {
    return RegisterChangeListenerForDataFilter(FindDataFilter (name), listener);
  }
  
  bool    UnregisterChangeListenerForDataFilter (_String const& name, BaseRef listener) {
    return UnregisterChangeListenerForDataFilter(FindDataFilter (name), listener);
  }
  
  long    FindDataFilter (_String const& name) {
    return _data_filters.Find (&name);
  }
  
  long    StoreDataFilter (_String const& name, _DataSetFilter* object, bool handle_errors) {
    
    if (name.IsValidIdentifier(fIDAllowCompound)) {
      long exists_already = FindDataFilter(name);
      
      /*printf ("[StoreDataFilter] %s %d\n", name.sData, exists_already);
       
       _SimpleList history;
       long locked_index = _data_filter_locks.Next (-1, history);
       while (locked_index >= 0) {
       printf ("\tLOCKED %s\n", GetFilterName((long)_data_filter_locks.Retrieve(locked_index))->sData);
       locked_index = _data_filter_locks.Next (locked_index, history);
       } */
      
      if (exists_already >= 0L) {
        if (_IsObjectLocked(exists_already, HY_BL_DATASET_FILTER)) {
          if (handle_errors) {
             hy_global::HandleApplicationError (_String ("DataSetFilter ") & name.Enquote() & " could not be created because an existing filter of the same name is locked");
          }
          return -1;
        }
        
        //DeleteObject ((_DataSetFilter*)_data_filters.GetXtra (exists_already));
        _data_filters.SetXtra(exists_already, object, false); // this will delete the existing object
        _NotifyDataFilterListeners (exists_already, kNotificationTypeChange);
        
      } else {
        exists_already = _data_filters.Insert (new _String(name), (long)object, false, false);
      }
      
      
      
      _SetDataFilterParameters (name, *object);
      return exists_already;
    } else {
      if (handle_errors) {
         hy_global::HandleApplicationError (_String ("The name ") & name.Enquote() & " is not a valid HyPhy id in call to " & __PRETTY_FUNCTION__);
      }
    }
    return -1;
  }
  
  bool    DeleteDataFilter (long index) {
    if (_data_filters.IsValidIndex (index)) {
      if (_IsObjectLocked(index, HY_BL_DATASET_FILTER)) {
        return false;
      }
      _NotifyDataFilterListeners (index, kNotificationTypeDelete);
      _KillDataFilterParameters( *GetFilterName (index));
      _data_filters.Delete ((BaseRef)GetFilterName(index), true);
    }
    return true;
  }
  
  bool    DeleteDataFilter (_String const& name) {
    return DeleteDataFilter (FindDataFilter(name));
  }
  
  
  _String const * GetFilterName (long index) {
    if (_data_filters.IsValidIndex(index)) {
      return (_String*)_data_filters.Retrieve(index);
    }
    return nil;
  }
  
  void ClearFilters (void) {
    // note that this ignores lock information
    _data_filters.Clear(true);
    _data_filter_locks.Clear ();
    
  }
  
  void ClearAllGlobals (void) {
    ClearFilters();
  }
  
  //____________________________________________________________________________________
  
  bool    IsModelReversible (long mid) {
    _Variable *m = nil,
              *f = nil;
    bool    mbf;
    RetrieveModelComponents (mid, m, f, mbf);
    if (m&&f) {
        return ((_Matrix*)m->GetValue())->IsReversible(mbf?nil:(_Matrix*)f->GetValue());
    }
    return false;
  }
  
  //____________________________________________________________________________________
  
  bool    IsModelOfExplicitForm (long modelID) {
    if (modelID != HY_NO_MODEL) {
      return modelTypeList.list_data[modelID] != 0;
    }
    return false;
  }
  
  //____________________________________________________________________________________
  
  _String const  GenerateUniqueObjectIDByType (_String const & base, const long type) {
    _AVLList   * names = nil;
    _List*       legacy_list = nil;
    
    switch (type) {
      case HY_BL_DATASET:
        legacy_list = &dataSetNamesList;
        break;
      case HY_BL_DATASET_FILTER:
        names = &_data_filters;
        break;
      case HY_BL_TREE:
        names = &variableNames;
        break;
    }
    
    _String try_name;
    
    if (names) {
      try_name = base;
      long suffix = 1L;
      
      while (names->Find (&try_name) >= 0) {
        try_name = base & "_" & suffix++;
      }
    } else if (legacy_list) {
      try_name = legacy_list->GenerateUniqueNameForList (base, false);
    }
    
    return try_name;
  }
  
  AVLListXLIterator  ObjectIndexer (const long type) {
    switch (type) {
      case HY_BL_DATASET_FILTER:
        return AVLListXLIterator (&_data_filters);
        break;
    }
    hy_global::HandleApplicationError (_String("Called ") & __PRETTY_FUNCTION__ & " with an unsupported type");
    return AVLListXLIterator (nil);
  }
  
  unsigned long  CountObjectsByType (const long type) {
    switch (type) {
      case HY_BL_DATASET_FILTER:
        return _data_filters.countitems();
        break;
    }
    hy_global::HandleApplicationError (_String("Called ") & __PRETTY_FUNCTION__ & " with an unsupported type");
    return 0UL;
  }
  
  
  //____________________________________________________________________________________
  
  _String const * GetObjectNameByType (const long type, const long index, bool correct_for_empties) {
    
    if (index < 0L) {
      return nil;
    }
    
    _List * theList = nil;
    switch (type) {
      case HY_BL_DATASET:
        theList = &dataSetNamesList;
        break;
      case HY_BL_DATASET_FILTER:
        theList = &_data_filter_aux;
        break;
        
      case HY_BL_LIKELIHOOD_FUNCTION:
        theList = &likeFuncNamesList;
        break;
      case HY_BL_HBL_FUNCTION:
        theList = &batchLanguageFunctionNames;
        break;
      case HY_BL_MODEL:
        theList = &modelNames;
        break;
      case HY_BL_SCFG:
        theList = &scfgNamesList;
        break;
      case HY_BL_BGM:
        theList = &bgmNamesList;
        break;
        
    }
    if (theList) {
      // account for deleted objects
      if (!correct_for_empties)
        return (_String*)theList->GetItemRangeCheck (index);
      
      long counter = 0;
      for (unsigned long name_index = 0; name_index < theList->lLength; name_index++) {
        _String *thisName = (_String*)theList->GetItem(name_index);
        if (thisName && !thisName->empty()) {
          if (name_index - counter == index) {
            return thisName;
          }
        } else {
          counter ++;
        }
      }
    }
    return nil;
  }
    
    //
    
    // SLATED for DEPRECATION
    
    //____________________________________________________________________________________
    
    //____________________________________________________________________________________
    
    BaseRefConst _HYRetrieveBLObjectByName    (_String const& name, long& type, long *index, bool errMsg, bool expression_based_lookup, _ExecutionList* current_program) {
        
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
            loc = FindDataFilter (name);
            if (loc >= 0) {
                type = HY_BL_DATASET_FILTER;
                if (index) {
                    *index = loc;
                }
                return GetDataFilter (loc);
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
            if (loc < 0 && name == hy_env::last_model_parameter_list || name == hy_env::use_last_model) {
                loc = lastMatrixDeclared;
            }
            if (loc >= 0) {
                type = HY_BL_MODEL;
                if (index) {
                    *index = loc;
                }
                if (IsModelOfExplicitForm(loc)) {
                    return (BaseRef)modelMatrixIndices.list_data[loc];
                }
                return LocateVar (modelMatrixIndices.list_data[loc]);
            }
        }
        
        if (type & HY_BL_HBL_FUNCTION) {
            loc = FindBFFunctionName(name);
            if (loc >= 0) {
                type = HY_BL_HBL_FUNCTION;
                if (index) {
                    *index = loc;
                }
                return &GetBFFunctionBody (loc);
            }
        }
        
        if (type & HY_BL_TREE) {
            _Variable* tree_var = FetchVar (LocateVarByName(name));
            if (tree_var && tree_var->ObjectClass() == TREE) {
                type = HY_BL_TREE;
                return tree_var;
            }
        }
        
        if (expression_based_lookup) {
            /**
             this is meant to resolve reference based expressions, like ^(string) or *(string)
             
             */
            
            _String referenced_object;
            
            hy_reference_type ref = name.ProcessVariableReferenceCases (referenced_object, (current_program && current_program->nameSpacePrefix) ? current_program->nameSpacePrefix->GetName() : nil);
            
            if (ref == kStringLocalDeference || ref == kStringGlobalDeference) {
                //printf ("%s\n", referenced_object.get_str());
                return _HYRetrieveBLObjectByName (referenced_object, type, index, errMsg, false, current_program);
            }
            
        }

        
        if (errMsg) {
            HandleApplicationError (_String ("'") & name & "' does not refer to an existing object of type " & _HYHBLTypeToText (type));
        }
        type = HY_BL_NOT_DEFINED;
        return nil;
    }
    
    //____________________________________________________________________________________
    

    BaseRef _HYRetrieveBLObjectByNameMutable    (_String const& name, long& type, long *index, bool errMsg, bool expression_based_lookup, _ExecutionList* current_program) {
        using namespace hyphy_global_objects;
        
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
            loc = FindDataFilter (name);
            if (loc >= 0) {
                type = HY_BL_DATASET_FILTER;
                if (index) {
                    *index = loc;
                }
                return ExclusiveLockDataFilter (loc);
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
            if (loc < 0 && (name == hy_env::last_model_parameter_list || name == hy_env::use_last_model)) {
                loc = lastMatrixDeclared;
            }
            if (loc >= 0) {
                type = HY_BL_MODEL;
                if (index) {
                    *index = loc;
                }
                if (IsModelOfExplicitForm(loc)) {
                    return (BaseRef)modelMatrixIndices.list_data[loc];
                }
                return LocateVar (modelMatrixIndices.list_data[loc]);
            }
        }
        
        if (type & HY_BL_HBL_FUNCTION) {
            loc = FindBFFunctionName(name);
            if (loc >= 0) {
                type = HY_BL_HBL_FUNCTION;
                if (index) {
                    *index = loc;
                }
                return &GetBFFunctionBody (loc);
            }
        }
        
        if (type & HY_BL_TREE) {
            _Variable* tree_var = FetchVar (LocateVarByName(name));
            if (tree_var && tree_var->ObjectClass() == TREE) {
                return tree_var;
            }
        }
        
        if (expression_based_lookup) {
            /**
             this is meant to resolve reference based expressions, like ^(string) or *(string)
             
             */
            
            _String referenced_object;
            
            hy_reference_type ref = name.ProcessVariableReferenceCases (referenced_object, (current_program && current_program->nameSpacePrefix) ? current_program->nameSpacePrefix->GetName() : nil);
            
            if (ref == kStringLocalDeference || ref == kStringGlobalDeference) {
                //printf ("%s\n", referenced_object.get_str());
                return _HYRetrieveBLObjectByNameMutable (referenced_object, type, index, errMsg, false, current_program);
            }
            
        }
        
        if (errMsg) {
            HandleApplicationError (_String ("'") & name & "' does not refer to an existing object of type " & _HYHBLTypeToText (type));
        }
        type = HY_BL_NOT_DEFINED;
        return nil;
    }
    
    //____________________________________________________________________________________
    long    FindDataSetName (_String const&s) {
        return dataSetNamesList.FindObject (&s);
    }
    
    //____________________________________________________________________________________
    long    FindBgmName (_String const&s) {
      return bgmNamesList.FindObject (&s);
    }

    //____________________________________________________________________________________
    long    FindScfgName (_String const&s) {
      return scfgNamesList.FindObject (&s);
    }
  
 
    //____________________________________________________________________________________
    long    FindBFFunctionName (_String const&s, _VariableContainer const* theP) {
      if (theP) {
        _String prefix = *(theP->GetName());
        
          //ReportWarning (_String ("Looking for ") & s.Enquote() & " in " & prefix.Enquote());
        
        while (1) {
          _String test_id = prefix & '.' & s;
          long idx = batchLanguageFunctionNamesIndexed.Find (&test_id);
          if (idx >= 0) {
            return batchLanguageFunctionNamesIndexed.GetXtra(idx);
              //s = test_id;
              //return idx;
          }
          long cut_at = prefix.FindBackwards ('.', 0, -1);
          if (cut_at > 0) {
            prefix.Trim (0, cut_at - 1);
          } else {
            break;
          }
        };
      }
      
      
        //ReportWarning (_String ("Looking for ") & s.Enquote() & " in global context");
      return batchLanguageFunctionNamesIndexed.FindAndGetXtra(&s,-1);
    }
    //____________________________________________________________________________________
    long    FindLikeFuncName (_String const&s, bool tryAsAString)
    {
        long try1 = likeFuncNamesList.FindObject (&s);
        if (try1 < 0 && tryAsAString) {
            _String s2 (ProcessLiteralArgument(&s, nil));
            try1 = likeFuncNamesList.FindObject(&s2);
        }
        return try1;
    }
    
    //____________________________________________________________________________________
    long    FindModelName (_String const &s) {
        if (s == hy_env::use_last_model) {
            return lastMatrixDeclared;
        }
        
        return modelNames.FindObject (&s);
    }
    
    //____________________________________________________________________________________
    _LikelihoodFunction*    FindLikeFuncByName (_String const&s)
    {
        long i = FindLikeFuncName(s);
        if (i>=0) {
            return (_LikelihoodFunction*)likeFuncList (i);
        }
        return nil;
    }
    
    //____________________________________________________________________________________
    long    FindLikeFuncIndex (void * const lfp) {
        return likeFuncList._SimpleList::Find ((long)lfp);
    }

    //____________________________________________________________________________________
    long    FindSCFGName (_String const&s)
    {
        return scfgNamesList.FindObject (&s);
    }

  
}
