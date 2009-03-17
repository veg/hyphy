/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006  
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

#ifndef _HLIST_
#define _HLIST_
//#pragma once
#include "baseobj.h"
#define  MEMORYSTEP 8

// for storing longs

//_________________________________________________________________________________________

class _SimpleList:public BaseObj{

			// contructor/destructor methods
			public:
			
			_SimpleList 					(); 
				//does nothing
			_SimpleList 					(unsigned long);
				//length constructor
			_SimpleList 					(_SimpleList&, long = 0, long = -1);
				// stack copy contructor
			_SimpleList 					(long);
				// data constructor (1 member list)
				
			_SimpleList						(long, long, long);
				// arithmetic series populator: size, first item, step
					
virtual 	~_SimpleList					(void);
 				//destructor

		 	long& operator []  				(long);
  				// element location functions - read/write
 			
			long operator () 				(unsigned long);
  				// element location functions - read only
  	  				
bool		Equal 							(_SimpleList&);

void		TrimMemory						(void);

void		RequestSpace 					(long);
				// request space for a given # of elements

virtual 	_SimpleList operator = 			(_SimpleList);
 				// assignment operator
 			
 			unsigned long 	countitems		(void);
 				// list length
 			
virtual 	_SimpleList operator & 			(_SimpleList);
 				// append operator
 			
virtual 	void operator << 				(long);
 				// append number to this
 				
virtual 	bool operator >>				(long);
 				// append number to this if it's not in the list (search first). List assumed unsorted.

virtual 	void operator << 				(_SimpleList&);
 			
virtual		void InsertElement 				(BaseRef br, long insertAt = -1, bool store = true, bool pointer = true);

  			void Clear	 					(void);
	
			long Element					(long);
			// much like [] and () except negative indices return offsets from the end
			// invalid indices return 0
 				 				
virtual 	long  Find 						(long, long startAt = 0);
 			// find the position of a search string in the list of strings (ONLY)
 			// -1 if not found
virtual		void  FilterRange				(long, long);
			// retain all those elements that are between (strictly) the 1st and the 2nd argument

virtual 	long  FindStepping				(long, long, long = 0);
 			
virtual 	long  BinaryFind 				(long, long startAt = 0);
 			// find the position of a search string in the list of strings (ONLY)
 			// -1 if not found

			long  BinaryInsert 		  		(long);
			// insert an element into the sorted list preserving the sortedness
			void  Delete 					(long, bool = true);
 			// delete the item at a given poisiton
 			
 			void  DeleteDuplicates 			(void);
 			// delete all duplicates in a sorted list

virtual		void  DeleteList 				(const _SimpleList&);

			void  Displace 					(long,long,long);
 			// shift a range of elements in the array

			void  Offset 					(long);
 			// add a number to each entry in the array

			bool  NChooseKInit 				(_SimpleList&, _SimpleList&, unsigned long, bool = false);
			// initialize the function to select all k-element subsets of a given simple list
			// arguments are:
			// [SimpleList] a state-storing simple list; will be approximately the same length as (*this) _SimpleList; 
			// DO NOT MANIPULATE this list outside NChooseKInit; it must persist between calls to NChooseK
			// [SimpleList] the receptacle list that will store k-tuples
			// [unsigned long]: how many elements to choose; must be <= lLength
			// [bool]: which algorithm to use for k-tuple generation; false - lexicographic (in the sense of the original list order)
			//		 : true - 'revolving door' method - TBA
			// return [bool]: true if successfully initialized
			
			bool NChooseK					(_SimpleList&, _SimpleList&);
			// select the next k-tuple
			// [SimpleList] the state-storing list previously populated by NChooseKInit
			// [SimpleList] the receptacle that will store k-tuples
			// returns [bool] true is more k-tuples are available; [false] if the last one has just been stored

			void  Permute 					(long);
			// permute elements in blocks of given size
			
			void  PermuteWithReplacement	(long);
			// permute elements in blocks of given size with possible replacement
			
			void  Union	   					(_SimpleList&, _SimpleList&);
			void  Intersect					(_SimpleList&, _SimpleList&);
			void  XOR						(_SimpleList&, _SimpleList&);
			void  Subtract					(_SimpleList&, _SimpleList&);
			long  CountCommonElements		(_SimpleList&, bool = false);
			 

			void  Merge 					(_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults = nil, _SimpleList* mergeResults2 = nil);

			void  Swap 						(long, long); //swap two elements 
			void  Populate					(long, long, long); // total elements,
																// start at
																// increment by

			void  Flip 						(void); //flip the order of list elements
	
virtual 	BaseRef toStr 					(void);

void		RecursiveIndexSort 				(long from, long to, _SimpleList* index);

BaseRef		ListToPartitionString 			(void);
 				
virtual		BaseRef	makeDynamic 			(void);

virtual		void	Initialize 				(bool = true);

virtual		void	Duplicate 				(BaseRef);

void				Sort 					(bool ascending = true);
void				BubbleSort 				(void);
void				QuickSort 				(long, long);
 				
long*		quickArrayAccess 				(void) 
												{ return (long*)lData;}
 				 			
virtual		long	Compare					(long,long);
virtual		long	Compare					(BaseRef,long);
	
static		void	NormalizeCoordinates    (long&, long&, const unsigned long);
			/* SLKP: 20090316
			 
				
			   given a range [from,to] and a given list,
			   make the range conform to the list (e.g. resolve negative to and/or from coordinates)
			   clip the range to fit the list etc
			 
               the third argument is the length of the list to normalize with respect to			 
			*/

			friend class _AVLList;

			protected:
 			
 			// data fields
 			
 			unsigned long laLength; 
 						  //memory allocated enough for this many slots

			public:
			
 			long* 		  lData;
  			unsigned long lLength;//actual length
			
};

//_________________________________________________________________________________________

class _List:public _SimpleList {
	
			// contructor/destructor methods
			public:
			
			_List (); 
				//does nothing
			_List (unsigned long);
				//length constructor
			_List (const _List&, long = 0, long = -1);
				// stack copy contructor
			_List (BaseRef, char);
				// construct a list of substrings from the original string separated by char
			_List (BaseRef);
				// data constructor (1 member list)
virtual 	~_List(void);
 				//destructor

			BaseRef& operator [] (long);
  				// element location functions - read/write
			
			BaseRef operator () (unsigned long);
  				// element location functions - read only

virtual 	_List operator = (_List&);
 				// assignment operator
 			 			
 			_List operator & (_List&);
 				// append operator
 			
 			void operator && (BaseRef);
 				// append duplicate to *this
 			
  			void operator && (const char*);
 				// append a string created from a literal

			void operator << (BaseRef);
 				// append reference to *this

			void operator << (_List&);
 				// append reference to *this

 			void AppendNewInstance (BaseRef);
 				// append reference to *this
 				
 			bool Equal		 (_List&);

			void Place (BaseRef);
 				// append reference to *this

			void  Replace (long, BaseRef,bool dup = true);
 			// replace an item 
		
virtual 	long	 FreeUpMemory (long);


 			void bumpNInst (void);
 				
virtual 	void Clear	   (void);

 			_List operator & (BaseRef);
 				// append operator
 				 				
virtual 	long  Find (BaseRef, long startat = 0);
 			// find the position of a search string in the list of strings (ONLY)
 			// -1 if not found

virtual 	long  FindString 		  (BaseRef, long startat = 0, bool caseSensitive = true, long upTo = -1);
 			// find the position of a search string in the list of strings (ONLY)
 			// -1 if not found
 			// faster than the one above, since it assumes string entries

virtual 	long  BinaryFind 		  (BaseRef );
 			// find the position of a search string in the list of strings (ONLY)
 			// uses binary sort on sorted lists
 			// -1 if not found

			long  BinaryInsert 		  (BaseRef);
			// insert an element into the sorted list preserving the sortedness

 			void  Delete (long );
 			// delete the item at a given poisiton

virtual		void  DeleteList (const _SimpleList&);
 			// delete the item at a given poisiton

virtual		void InsertElement (BaseRef br, long insertAt = -1, bool store = true);
 

virtual 	BaseRef toStr 			(void);
 				
virtual		void	toFileStr 		(FILE*);

virtual		BaseRef	makeDynamic		(void);

virtual		void	Duplicate 		(const BaseRef);

virtual		long	Compare			(long,long);
virtual		long	Compare			(BaseRef,long);

			void  	Intersect		(_List&, _List&, _SimpleList* = nil, _SimpleList* = nil);
				
 				 			 			
};

//_________________________________________________________________________________________

class _AVLList: public BaseObj {

public:

			_AVLList	(_SimpleList*);
virtual		~_AVLList	(void) {}


	long				Find 				(BaseRef);
	char				FindBest			(BaseRef, long&);
	long				Find 				(BaseRef, _SimpleList&);
	long				Next 				(long,    _SimpleList&);
	long				Prev 				(long,    _SimpleList&);
	long				First				(void);
	long				Last				(void);
	long				Insert				(BaseRef, long = 0, bool = true);
											// the bool flag is to say whether to dup the object being inserted
	void				Delete				(BaseRef, bool = false);
	virtual	 void		ReorderList			(_SimpleList* = nil);
	BaseRef				Retrieve			(long);
	virtual	 long	 	InsertData			(BaseRef, long, bool);
	unsigned long		countitems			(void);
	virtual  void		Clear				(bool = false);
	virtual  BaseRef	toStr				(void);
	virtual  bool		HasData				(long);
	virtual  long		Traverser			(_SimpleList&, long &, long = -1);
	virtual  long		GetRoot				(void) { return root;}
	virtual  void		DeleteXtra			(long) {};
	virtual  void		DeleteAll			(bool cL) {Clear(cL); DeleteObject (dataList); }
	
	// data members
	
	void				ConsistencyCheck    (void);
	
	_SimpleList *dataList,
				leftChild,
				rightChild,
				balanceFactor,
				emptySlots;
				
	long		root;

};

//_________________________________________________________________________________________

class _AVLListX: public _AVLList {

public:

			_AVLListX	(_SimpleList*);
virtual		~_AVLListX	(void) {}


	
	virtual  BaseRef	toStr				(void);
	virtual	 long	 	InsertData			(BaseRef, long, bool);
			 long		GetXtra				(long);
			 void		SetXtra				(long,long);
	virtual  void		Clear				(bool = false);
	virtual  void		DeleteXtra			(long);
				 
	_SimpleList			xtraD;
	
};

//_________________________________________________________________________________________

class _AVLListXL: public _AVLList {

public:

			_AVLListXL	(_SimpleList*);
virtual		~_AVLListXL	(void) {}


	virtual	 long	 	InsertData			(BaseRef, long,bool);
			 BaseRef	GetXtra				(long);
			 void		SetXtra				(long,BaseRef,bool);
	virtual  void		Clear				(bool = false);
	virtual  void		DeleteXtra			(long);
			 
	_List				xtraD;
	
};

//_________________________________________________________________________________________
 			
void	    SortLists (_SimpleList*, _SimpleList*);

extern		_List	  pathNames;
 			
#endif
 			