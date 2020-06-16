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


#ifndef __GLOBAL_OBJECT_LISTS_H__

#define  __GLOBAL_OBJECT_LISTS_H__
#include "dataset_filter.h"
#include "avllistxl_iterator.h"

/**
 The namespace and attendant functions are meant as an intermediate
 phase for controlling access to non-variable globally accessible types:
 DataSet, DataSetFilter, LikelihoodFunction, SCFG, BGM, Model
 */

namespace hyphy_global_objects {
  
  /**
   Functions for accessing and manipulating DataSetFilters
   */
  
  const   _DataSetFilter * GetDataFilter (long index);
  /**
   Retrieve a read-only pointer to a DataSetFilter object
   
   @param index the direct referencing index for this object
   @return the object pointer or null if reference is invalid
   
   @see get_data_filter (_String const)
   */
  
  const   _DataSetFilter * GetDataFilter (_String const& name);
  
  /**
   Retrieve a read-only pointer to a DataSetFilter object
   
   @name the globally scoped name of this object
   @return the object pointer or null if the name is invalid
   
   @see get_data_filter (long index)
   */
  
  _DataSetFilter * ExclusiveLockDataFilter (long index);
  
  /**
   Retrieve a read-write reference to a DataSetFilter object and mark it
   as exclusively locked; trying to delete this object or exclusively
   lock it again will return a null
   
   @param index the direct referencing index for this object
   @return the object pointer or null if reference is invalid or the lock could\
   not be obtained
   
   @see get_data_filter (_String const)
   @see release_data_filter_lock
   */
  
  _DataSetFilter * ExclusiveLockDataFilter (_String const& name);
  
  /**
   Retrieve a read-write reference to a DataSetFilter object and mark it
   as exclusively locked; trying to delete this object or exclusively
   lock it again will return a null
   
   @name the globally scoped name of this object
   @return the object pointer or null if reference is invalid or the lock could\
   not be obtained
   
   @see get_data_filter (_String const)
   @see release_data_filter_lock
   */
  
  bool    ReleaseDataFilterLock (long index);
  /**
   Release a previously acquired lock to this data set filter object
   
   @param index the direct referencing index for this object
   @return true if the object is now unlocked
   
   @see exclusive_lock_data_filter (_String const)
   @see release_data_filter_lock
   */
  bool    ReleaseDataFilterLock (_String const& name);
  /**
   Release a previously acquired lock to this data set filter object
   
   @param name the globally scoped name of this object
   @return true if the object is now unlocked
   
   @see exclusive_lock_data_filter (_String const)
   @see release_data_filter_lock
   */
  
  long    StoreDataFilter  (_String const& name, _DataSetFilter *filter, bool handle_errors = true);
  
  /**
   Store a reference to a DataSetFilter object in the global namespace
   
   @paramename the globally scoped name of this object
   @param filter the object to store; this operation will NOT increase a reference count of filter
   @param handle_errors if true, calls error handling functions if something bad happens
   @return the index of the stored object or <0 if the store operation failed
   
   */
  
  const _String* GetFilterName  (long index);
  
  /**
   Retrieve the name of the filter associated with the index
   
   @param index the direct referencing index for this object
   @return the name of the associated object, or NULL if failed (no such index)
   
   */
  
  long    FindDataFilter   (_String const& name);
  
  /**
   Find the referencing index for a given DataSetFilter object name
   @name the globally scoped name of this object
   @return the index of the corresponding object or <0 if no such object exists
   
   */
  
  bool    DeleteDataFilter (long index);
  
  /**
   Delete the DataSetFilter object
   
   @param index the direct referencing index for this object
   @return true if the object was deleted
   
   */
  
  bool    DeleteDataFilter (_String const& name);
  
  /**
   Delete the DataSetFilter object
   
   @param name the globally scoped name of this object
   @return true if the object was deleted
   
   */
  
  void    ClearFilters (void);
  
  /**
   Delete all DataSetFilter objects
   
   */
  
  void    ClearAllGlobals (void);
  /**
   Delete all global objects
   
   */
  
  
  bool    RegisterChangeListenerForDataFilter (const long index, BaseRef listener);
  /**
   Add an object to the list of listeners that will be notified if a data filter changes
   
   @param index the index of the datafilter to which the listener will be bound
   @param listener the object that will be listening to changes on this filter
   
   @return whether or not the registration was successful
   
   */
  
  bool    RegisterChangeListenerForDataFilter (_String const & name, BaseRef listener);
  /**
   Add an object to the list of listeners that will be notified if a data filter changes
   
   @param the name of the datafilter to which the listener will be bound
   @param listener the object that will be listening to changes on this filter
   
   @return whether or not the registration was successful
   
   */
  

  bool    UnregisterChangeListenerForDataFilter (const long index, BaseRef listener);
  /**
   Remove an object from the list of listeners that will be notified if a data filter changes
   
   @param index the index of the datafilter from which the listener will be unbound
   @param listener the object that will no longer be listening to changes on this filter
   
   @return whether or not the un-registration was successful
   
   */
  
  bool    UnregisterChangeListenerForDataFilter (_String const & name, BaseRef listener);
  /**
   Remove an object from the list of listeners that will be notified if a data filter changes
   
   @param the name of the datafilter from which the listener will be unbound
   @param listener the object that will no longer be listening to changes on this filter
   
   @return whether or not the un-registration was successful
   
   */

  unsigned long CountObjectsByType (const long type);
  /**
   Returns the count of objects for a specific type
   
   @param type the type of the object
   
   @return the number of objects of the requested type
   
   */
  
  
  AVLListXLIterator  ObjectIndexer (const long type);
  /**
   Returns a range-compiant object to iterate over the requested object type
   
   @param type the type of the object
   
   @return an iterator for the right object type
   
   */
  
  
  const   _String* GetObjectNameByType       (const long type, const long index, bool correct_for_empties = true);
  
  bool    IsModelReversible                 (long model_index);
  /** is the model reversible?
    @param model_index : model index
   
    @return T/F
   */
  bool    IsModelOfExplicitForm             (long model_index);
  /** is the model specified in a form that doesn't require exponentiation
   @param model_index : model index
   @return T/F
   */

  
  const   _String  GenerateUniqueObjectIDByType (_String const & base, const long type);
  /**
   Generate a unique object ID of a given type of the form base_[autogenerated number]
   
   @param base the prefix to use for the name
   @type the type of the object
   
   @return a unique ID
   
   */
    
    // TODO: 20171005 SLKP, these will be deprecated
    
    BaseRefConst _HYRetrieveBLObjectByName              (_String const& name, long& type, long* index = nil, bool errMsg = false, bool tryLiteralLookup = false, _ExecutionList* current_program = nil);
    
    BaseRef      _HYRetrieveBLObjectByNameMutable       (_String const& name, long& type, long* index = nil, bool errMsg = false, bool tryLiteralLookup = false, _ExecutionList* current_program = nil);
    

    long      FindDataSetName                 (_String const&);
    long      FindSCFGName                    (_String const&);
    long      FindBFFunctionName              (_String const&, _VariableContainer const* = nil);
    long      FindBgmName                     (_String const&);
    // added by afyp, March 18, 2007

    long      FindLikeFuncName                (_String const&, bool = false);
    long      FindModelName                   (_String const&);

    extern   _AVLListX batchLanguageFunctionNamesIndexed;
  
}

#endif
