/*
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
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


/** Legacy declarations */

extern _List batchLanguageFunctionNames;


_List        _data_filter_aux;

_AVLListXL   _data_filters (&_data_filter_aux);

_SimpleList  _data_filter_locks_aux;

_AVLList     _data_filter_locks (&_data_filter_locks_aux);



namespace hyphy_global_objects {
  
  
  
  /** 
      data filter objects -- internal
   */
  
  
  
  void  _SetDataFilterParameters (_String const& name, _DataSetFilter const& filter) {
    setParameter (WrapInNamespace ("species", &name), filter.NumberSpecies());
    setParameter (WrapInNamespace ("sites", &name), filter.GetSiteCount()/filter.GetUnitLength());
    
    _Parameter size_cutoff;
    checkParameter  (defaultLargeFileCutoff,size_cutoff, 100000.);
    
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
      return source->FindLong (index);
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

  
  long    FindDataFilter (_String const& name) {
    return _data_filters.Find (&name);
  }

  long    StoreDataFilter (_String const& name, _DataSetFilter* object, bool handle_errors) {
    
    if (name.IsValidIdentifier(true)) {
      long exists_already = FindDataFilter(name);
      if (exists_already >= 0L) {
        if (_IsObjectLocked(exists_already, HY_BL_DATASET_FILTER)) {
          if (handle_errors) {
            WarnError (_String ("DataSetFilter ") & name.Enquote() & " could not be stored because the existing filter of the same name is locked");
          }
          return -1;
        }
        DeleteObject ((_DataSetFilter*)_data_filters.GetXtra (exists_already));
        _data_filters.SetXtra(exists_already, object, false);
        
      } else {
        exists_already = _data_filters.Insert (new _String(name), (long)object, false, false);
      }
      

      _SetDataFilterParameters (name, *object);
      return exists_already;
    } else {
      if (handle_errors) {
        WarnError (_String ("The name ") & name.Enquote() & " is not a valid HyPhy id in call to (store_data_filter)");
      }
    }
    return -1;
  }
  
  bool    DeleteDataFilter (long index) {
    if (_data_filters.IsValidIndex (index)) {
      if (_IsObjectLocked(index, HY_BL_DATASET_FILTER)) {
        return false;
      }
      _KillDataFilterParameters( *GetFilterName (index));
     _data_filters.Delete ((BaseRef)index, true);
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
    WarnError (_String("Called ") & __PRETTY_FUNCTION__ & " with an unsupported type");
    return AVLListXLIterator (nil);
  }

  unsigned long  CountObjectsByType (const long type) {
    switch (type) {
      case HY_BL_DATASET_FILTER:
        return _data_filters.countitems();
        break;
    }
    WarnError (_String("Called ") & __PRETTY_FUNCTION__ & " with an unsupported type");
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
        if (thisName && thisName->sLength) {
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
  
}
