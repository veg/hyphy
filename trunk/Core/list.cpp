/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    		   (apoon@cfenet.ubc.ca)

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

_AVLList	structure inspired by the excellent documentation of 
GNU libavl 2.0.1 by Ben Pfaff (http://www.msu.edu/~pfaffben/avl/index.html)

*/

#include "hy_strings.h"
#include "errorfns.h"
#include "hy_lists.h"
#include "parser.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif


//______________________________________________________________
//does nothing

_SimpleList::_SimpleList ()
{
	Initialize(false);
	
}

//______________________________________________________________
//does nothing

_List::_List ()
{
}

//______________________________________________________________

void	_SimpleList::Initialize (bool doMemAlloc)
{
	BaseObj::Initialize();
	lLength = 0;
	if (doMemAlloc)
	{
		laLength = MEMORYSTEP;
		lData = (long*)MemAllocate (laLength * sizeof(Ptr));
	}
	else
	{
		laLength = 0;
		lData    = nil;
	}
}

//______________________________________________________________

void	_SimpleList::Duplicate (BaseRef theRef)
{
	_SimpleList* l	= (_SimpleList*)theRef;
	lLength			= l->lLength;
	laLength		= l->laLength;
	lData			= l->lData;
	if (lData) {
		checkPointer (lData = (long*)MemAllocate (laLength*sizeof (Ptr)));
		memcpy ((char*)lData, (char*)l->lData, lLength*sizeof (Ptr));
	}
}

//______________________________________________________________

void	_List::Duplicate (const BaseRef theRef)
{
	_SimpleList::Duplicate (theRef);
	if (lData)
		for (unsigned long i = 0; i<lLength; i++)
		{
			if (((BaseRef*)lData)[i])
				(((BaseRef*)lData)[i])->nInstances++;
		}
}
		
//______________________________________________________________
//length constructor

_SimpleList::_SimpleList (unsigned long l)
{
	lLength = 0;
	laLength = (l/MEMORYSTEP + 1)*MEMORYSTEP;
	lData = (long*)MemAllocate (laLength * sizeof(Ptr));
	memset (lData,0,laLength * sizeof(Ptr));

}

//______________________________________________________________
//length constructor and populator

_SimpleList::_SimpleList (long l, long start, long step)
{
	Initialize (false);
	Populate   (l,start,step);
}

//______________________________________________________________
//length constructor and populator

void _SimpleList::Populate (long l, long start, long step)
{
	RequestSpace (l);
	for (long k = 0; k < l; k++, start+=step)
		lData[k] = start;

	lLength = l;
}

//______________________________________________________________
//length constructor

_List::_List (unsigned long l):_SimpleList(l)
{
}
				
//______________________________________________________________
// stack copy contructor
_List::_List (const _List& l, long from, long to)
{
	if (from == 0 && to == -1) // copy the whole thing
	{
		BaseRef	  br = (BaseRef)&l;
		Duplicate (br);
	}
	else
	{
		Initialize			 ();
		NormalizeCoordinates (from, to, l.lLength);
			
		for (long i = from; i <= to; i++)
			(*this) << ((BaseRef*)l.lData)[i];	
	}	
}

//______________________________________________________________
// stack copy contructor
_List::_List (BaseRef ss, char sep)
{
	_String* s = (_String*)ss;
	if (s->Length()!=0) 
	{
		long cp=0,cpp;
		while ((cpp = s->Find(sep,cp,-1))!=-1)
		{
			if (cpp>cp)
				AppendNewInstance (new _String(*s,cp,cpp-1));
			else
				AppendNewInstance (new _String);
			cp=cpp+1;
		}
		
		AppendNewInstance (new _String(*s,cp,-1));
	}
}


//______________________________________________________________
// coordinate normalizer 
void _SimpleList::NormalizeCoordinates (long& from, long& to, const unsigned long refLength)
{
	if (to < 0)
		to = refLength+to;
	else
		to = to < refLength-1 ? to : refLength - 1;
	if (from < 0)
		from = refLength+from;
}

//______________________________________________________________
// stack copy contructor
_SimpleList::_SimpleList (_SimpleList& l, long from, long to)
{
	if (from == 0 && to == -1) // copy the whole thing
		Duplicate (&l);
	else
	{
		Initialize			 ();
		NormalizeCoordinates (from, to, l.lLength);
		
		for (long i = from; i < to; i++)
			(*this) << l.lData[i];	
	}
}

//______________________________________________________________
void	_SimpleList::Offset (long shift)
{
	if (lData)
	{
		for (long k=0; k<lLength; k++)
			lData[k] += shift;
	}
}

//______________________________________________________________
// data constructor (1 member list)
_List::_List (BaseRef br)
{
	lLength = 1;
	laLength = MEMORYSTEP;
	lData = (long*)MemAllocate (laLength * sizeof(Ptr));
	((BaseRef*)lData)[0]= br->makeDynamic();
}

//______________________________________________________________
// data constructor (1 member list)
_SimpleList::_SimpleList (long br)
{
	lLength = 1;
	laLength = MEMORYSTEP;
	lData = (long*)MemAllocate (laLength * sizeof(Ptr));
	((long*)lData)[0]= br;
}
	
//______________________________________________________________
_SimpleList::~_SimpleList(void)
//destructor
{
	if (nInstances<=1)
	{
		if (lData) 
			free (lData);
	}
	else
		nInstances--;	


}
//______________________________________________________________
_List::~_List(void)
//destructor
{
	if (nInstances<=1)
	{
		for (unsigned long i = 0; i<lLength; i++)
		{
			BaseRef t = ((BaseRef*)lData)[i];
			if (t){
				if (t->nInstances<=1){
					DeleteObject(t);
				}
				else t->nInstances--;
			}
		}
	}


}

//______________________________________________________________
// element location functions (0,llength - 1)
// negative indices return offsets from the end of the list

long _SimpleList::Element (long index)
{
	if (index >= 0)
	{
		if (index < lLength)
			return lData[index];
	}
	else
	{
		if (-index <= lLength)
			return lData[(long)lLength+index];
	}
	return 0;
}

//______________________________________________________________

long _SimpleList::Pop (void)
{
	if (lLength > 0)
	{
		lLength --;
		return lData[lLength];
	}
		
	return 0;
}



//______________________________________________________________
// element location functions (0,llength - 1)

long& _SimpleList::operator [] (long i)
{
	if (lLength == 0) 
		return lData[0];
	
	unsigned long in = (unsigned long)i;
	if (in>lLength-1) 
		in = lLength-1;
	
	return lData[in];
}

//______________________________________________________________
// element location functions (0,llength - 1)

BaseRef& _List::operator [] (long i)
{
	BaseRef t = BaseRef(_SimpleList::operator[] (i));
	if (t)
		if (t->nInstances>1)
		{
			t->nInstances--;
			((BaseRef*)(lData))[i]=t->makeDynamic();
		}
	
	return ((BaseRef*)(lData))[i];
}
//______________________________________________________________
// element location functions (0,llength - 1)

long _SimpleList::operator () (unsigned long i)
{
	//if (lLength == 0) return 0;
	//if (i>=lLength) i = lLength-1;
	return lData[i];
}

//______________________________________________________________
// element location functions (0,llength - 1)

BaseRef _List::operator () (unsigned long i)
{
	return ((BaseRef*)lData)[i];
}


//______________________________________________________________
// assignment operator
_SimpleList _SimpleList::operator = (_SimpleList l)
{
	Clear();
	lLength  = l.lLength;
	laLength = l.laLength;
	if (laLength)
	{
		checkPointer (lData = (long*)MemAllocate (laLength*sizeof (Ptr)));
		if (lLength)
			memcpy (lData,l.lData,lLength*sizeof (Ptr));
	}
	
	return *this;
}

//______________________________________________________________
// assignment operator
_List _List::operator = (_List& l)
{
	this->~_List();
	lLength = l.lLength;
	laLength = l.laLength;
	lData = l.lData;
	l.nInstances++;
	for (unsigned long i = 0; i<lLength; i++)
	{
		((BaseRef*)(lData))[i]->nInstances++;
	}
	return *this;
}

//______________________________________________________________
// list length
unsigned long _SimpleList::countitems(void)
{
	return lLength;
}

//______________________________________________________________

bool _SimpleList::Equal(_SimpleList& l2)
{
	if (lLength!=l2.lLength) 
		return false;
		
	for (long i=0;i<lLength; i++)
		if (lData[i] != l2.lData[i]) 
			return false;
			
	return true;
}

//______________________________________________________________

bool _List::Equal(_List& l2)
{
	if (lLength!=l2.lLength) 
		return false;
		
		
	for (long i=0;i<lLength; i++)
		if (!((_String*)lData[i])->Equal ((_String*)l2.lData[i]))
			return false; 
			
	return true;
}

//______________________________________________________________
// merge 2 lists (sorted)
void _SimpleList::Merge(_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults1, _SimpleList* mergeResults2)
{
	Clear();
	if (mergeResults1) 
		mergeResults1->Clear();
	if (mergeResults2) 
		mergeResults2->Clear();
	char    advancing = -1;
	long	* list1 = l1.quickArrayAccess(), 
			* list2 = l2.quickArrayAccess(),
			pos1=0, 
			pos2 =0, 
			nt1 = l1.lLength, 
			nt2 = l2.lLength,
			c,
			i;
	
	if (mergeResults1 && mergeResults2)
	{
		bool   doMerge1 = false, 
			   doMerge2 = false;
		while (1) // stuff left to do
		{
			if (advancing == 0) // advancing in the 1st list
			{
				pos1++;
				if (pos1==nt1)
				{
					advancing = 2;
					continue;
				}
				list1++;
				c = *list1-*list2;
				if (c<=0)
				{
					if (doMerge1)
						(*mergeResults1)<<lLength;
					
					(*this)<<*list1;
					if (c<0) 
					{
						if (!doMerge2)
						{
							if (pos1>=pos2)
							{
								doMerge2 = true;
								for (i=0; i<pos2; i++)
									(*mergeResults2)<<i;
							}
						}
						continue;
					}
				}
				if (c>0)
				{
					advancing = 1;
					if (!doMerge1)
					{
						for (i=0; i<pos1; i++)
							(*mergeResults1)<<i;
						doMerge1 = true;
					}
					if (doMerge2)
						(*mergeResults2)<<lLength;
					(*this)<<*list2;
					continue;
				}
					
			}
			else
				if (advancing == 1) // advancing in the 2nd list
				{
					pos2++;
					if (pos2==nt2)
					{
						advancing = 3;
						continue;
					}
					list2++;
					c = *list2-*list1;
					if (c<=0)
					{
						if (doMerge2)
							(*mergeResults2)<<lLength;
						(*this)<<*list2;
						if (c<0) 
						{
							if (!doMerge1)
							{
								if (pos2>=pos1)
								{
									doMerge1 = true;
									for (i=0; i<pos1; i++)
										(*mergeResults1)<<i;
								}
							}
							continue;
						}
					}
					if (c>0)
					{
						advancing = 0;
						if (!doMerge2)
						{
							for (i=0; i<pos2; i++)
								(*mergeResults2)<<i;
							doMerge2 = true;
						}
						if (doMerge1)
							(*mergeResults1)<<lLength;
						(*this)<<*list1;
						continue;
					}
				}
				else
					if (advancing == 2) // flush out the 2nd list
					{
						if (!doMerge1&&(pos2<nt2))
						{
							for (i=0; i<nt1; i++)
								(*mergeResults1)<<i;
						}
						if (doMerge2)
							while (pos2<nt2)
							{	
								(*mergeResults2)<<lLength;
								(*this)<<*list2;
								list2++;
								pos2++;
							}
						else
							while (pos2<nt2)
							{	
								(*this)<<*list2;
								list2++;
								pos2++;
							}
						break;
					}
					else
						if (advancing == 3) // flush out the 1st list
						{
							if (!doMerge2&&(pos1<nt1))
							{
								for (i=0; i<nt2; i++)
									(*mergeResults2)<<i;
							}
							if (doMerge1)
								while (pos1<nt1)
								{	
									(*mergeResults1)<<lLength;
									(*this)<<*list1;
									list1++;
									pos1++;
								}
							else
								while (pos1<nt1)
								{	
									(*this)<<*list1;
									list1++;
									pos1++;
								}
							break;
						}
						else
							if (advancing == -1) // just starting
							{
								if (!nt1) // first list is empty!
								{
									advancing = 2;
									continue;
								}
								if (!nt2) // second list is empty!
								{
									advancing = 3;
									continue;
								}
								c = *list1-*list2;
								if (c<=0) // begin with the first list
								{
									(*this)<<*list1;
									advancing = 0;
									if (c)
									{
										doMerge2 = true;
										continue;
									}
								}
								else
								{
									(*this)<<*list2;
									advancing = 1;
									doMerge1 = true;
									continue;
								}
									
							}
							
				if (advancing == 0) // moving up in the second term
				{
					pos1++;
					if (pos1==nt1)
					{
						list2++;
						pos2++;
						if (doMerge2)
						{
							(*mergeResults2)<<lLength-1;
						}
						advancing = 2;
						continue;
					}
					else
					{
						advancing = 1;
						if (doMerge2)
						{
							(*mergeResults2)<<lLength-1;
						}
						list1++;
					}
				}
				else
				{
					pos2++;
					if (pos2==nt2)
					{
						list1++;
						pos1++;
						if (doMerge1)
						{
							(*mergeResults1)<<lLength-1;
						}
						advancing = 3;
						continue;
					}
					else
					{
						list2++;
						if (doMerge1)
						{
							(*mergeResults1)<<lLength-1;
						}
						advancing = 0;
					}
				}
			}
	}
	else
	{
		while (1) // stuff left to do
		{
			if (advancing == 0) // advancing in the 1st list
			{
				pos1++;
				if (pos1==nt1)
				{
					advancing = 2;
					continue;
				}
				list1++;
				c = *list1-*list2;
				if (c<=0)
				{
					(*this)<<*list1;
					if (c<0) continue;
				}
				if (c>0)
				{
					advancing = 1;
					(*this)<<*list2;
					continue;
				}
					
			}
			else
				if (advancing == 1) // advancing in the 2nd list
				{
					pos2++;
					if (pos2==nt2)
					{
						advancing = 3;
						continue;
					}
					list2++;
					c = *list2-*list1;
					if (c<=0)
					{
						(*this)<<*list2;
						if (c<0) continue;
					}
					if (c>0)
					{
						advancing = 0;
						(*this)<<*list1;
						continue;
					}
				}
				else
					if (advancing == 2) // flush out the 2nd list
					{
						while (pos2<nt2)
						{	
							(*this)<<*list2;
							list2++;
							pos2++;
						}
						break;
					}
					else
						if (advancing == 3) // flush out the 2nd list
						{
							while (pos1<nt1)
							{	
								(*this)<<*list1;
								list1++;
								pos1++;
							}
							break;
						}
						else
							if (advancing == -1) // just starting
							{
								if (!nt1) // first list is empty!
								{
									advancing = 2;
									continue;
								}
								if (!nt2) // second list is empty!
								{
									advancing = 3;
									continue;
								}
								c = *list1-*list2;
								if (c<=0) // begin with the first list
								{
									(*this)<<*list1;
									advancing = 0;
									if (c)
									{
										continue;
									}
								}
								else
								{
									(*this)<<*list2;
									advancing = 1;
									continue;
								}
									
							}
							
				if (advancing == 0) // moving up in the second term
				{
					pos1++;
					if (pos1==nt1)
					{
						list2++;
						pos2++;
						advancing = 2;
						continue;
					}
					else
					{
						advancing = 1;
						list1++;
					}
				}
				else
				{
					pos2++;
					if (pos2==nt2)
					{
						list1++;
						pos1++;
						advancing = 3;
						continue;
					}
					else
					{
						list2++;
						advancing = 0;
					}
				}
				
		
		}
	}
}
 			
//______________________________________________________________
// append operator
_List _List::operator & (_List& l)
{
	_List res (l.lLength + lLength);
	if (!res.laLength) return res;
	
	if (lData&&lLength) {
		memcpy(res.lData,lData,lLength*sizeof(void*));
	}
	if (l.lData&&l.lLength) {
		memcpy((char*)res.lData+lLength*sizeof(void*),l.lData, l.lLength*sizeof(void*));
	}
	res.lLength = l.lLength + lLength;
	unsigned long i;
	for (i = 0; (i<lLength); i++)
	{
		((BaseRef*)lData)[i]->nInstances++;

	}
	for (i=0; i<l.lLength; i++)
	{
		((BaseRef*)l.lData)[i]->nInstances++;
	}
	return res;
}

//______________________________________________________________
// append operator
_SimpleList _SimpleList::operator & (_SimpleList l)
{
	_SimpleList res (l.lLength + lLength);
	if (!res.laLength) return res;
	
	if (lData&&lLength) {
		memcpy(res.lData,lData,lLength*sizeof(void*));
	}
	if (l.lData&&l.lLength) {
		memcpy((char*)res.lData+lLength*sizeof(void*),l.lData, l.lLength*sizeof(void*));
	}
	res.lLength = l.lLength + lLength;

	return res;
}

//______________________________________________________________
// append operator
_List _List::operator & (BaseRef br)
{
	_List res (lLength+1);
	if (!res.laLength) return res;
	
	if (lData) {
		memcpy(res.lData,lData,lLength*sizeof(void*));
	}
	for (unsigned long i = 0; (i<lLength); i++)
	{
		((BaseRef*)lData)[i]->nInstances++;
	}
	res.lLength=lLength+1;
	((BaseRef*)res.lData)[lLength]=br->makeDynamic();
	return res;
}

//______________________________________________________________
void _List::operator && (BaseRef br)
{
	InsertElement (br);
}

//______________________________________________________________
void _List::operator && (const char * buffer)
{
	_String* 		newString = new _String (buffer);
	checkPointer    (newString);
	InsertElement   (newString,-1,false);
	DeleteObject    (newString);
}


//______________________________________________________________
void _List::operator << (BaseRef br)
{
//	InsertElement (br, -1, false);
	lLength++;
	if (lLength>laLength)
	{
		unsigned long incBy = (MEMORYSTEP*5 > lLength)? MEMORYSTEP: lLength/5;
		
		laLength+=incBy;
		
		if (lData)
			checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
		else
			checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
	}
	((BaseRef*)lData)[lLength-1]=br;
	br->nInstances++;
}

//______________________________________________________________
void _List::AppendNewInstance (BaseRef br)
{
	if (br)
	{
		(*this)<<br;
		br->nInstances--;
	}
	else
		checkPointer (br);
}


//______________________________________________________________
void _List::operator << (_List& source)
{
	for (long k=0; k<source.lLength; k++)
		(*this) << ((BaseRef*)source.lData)[k];
}

//______________________________________________________________
void _List::Place (BaseRef br)
{
//	InsertElement (br, -1, false);
	lLength++;
	if (lLength>laLength)
	{
		laLength+=MEMORYSTEP;
		if (lData)
			checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
		else
			checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
	}
	((BaseRef*)lData)[lLength-1]=br;
}

//______________________________________________________________
void _SimpleList::operator << (long br)
{
	InsertElement ((BaseRef)br, -1, false, false);
}

//______________________________________________________________
bool _SimpleList::operator >> (long br)
{
	if (Find(br) == -1)
	{
		InsertElement ((BaseRef)br, -1, false, false);
		return true;
	}
	return false;
}

//______________________________________________________________
void _SimpleList::operator << (_SimpleList& source)
{
	for (long k=0; k<source.lLength; k++)
		(*this) << source.lData[k];
}

//______________________________________________________________
// append & store operator
void _SimpleList::InsertElement (BaseRef br, long insertAt, bool store, bool pointer)
{
	lLength++;
	if (lLength>laLength)
	{
		unsigned long incBy = (MEMORYSTEP*5 > lLength)? MEMORYSTEP: lLength/5;
		
		laLength+=incBy;

		//memAlloc += sizeof(Ptr)*incBy;

		if (lData)
			lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*));
		else
			lData = (long*)MemAllocate(laLength*sizeof(void*));
			
		if (!lData)
			checkPointer (lData);
	}
	if (insertAt==-1)
	{
		if (store)
			((BaseRef*)lData)[lLength-1]=br->makeDynamic();
		else
		{
			((BaseRef*)lData)[lLength-1]=br;
			if (pointer)
				br->nInstances++;
		}
	}
	else
	{
		//insertAt = insertAt>=lLength?lLength:insertAt;
		insertAt = insertAt>=lLength?lLength-1:insertAt;
		long     moveThisMany = (laLength-insertAt-1);
		if (moveThisMany < 32)
			for (long k=insertAt+moveThisMany; k> insertAt ; k--)
				lData[k] = lData[k-1];
		else
		{
			memmove (((char**)lData)+(insertAt+1), ((char**)lData)+insertAt, moveThisMany*sizeof(void*));
		}
		
		if (store)
			((BaseRef*)lData)[insertAt]=br->makeDynamic();
		else
		{
			((BaseRef*)lData)[insertAt]=br;
			if (pointer)
				br->nInstances++;
		}
	}
		

}

//______________________________________________________________
// convert a list into a partition style string
BaseRef _SimpleList::ListToPartitionString ()
{
	_String *result = new _String ((unsigned long)64,true),
			conv;
	
	for (long k=0; k<lLength; k++)
	{
		long m;
		for (m=k+1; m<lLength; m++)
			if (lData[m]-lData[m-1]!=1)
				break;
		if (m>k+2)
		{
			conv = lData[k];
			(*result) << & conv;
			(*result) << '-';
			conv = lData[m-1];
			(*result) << & conv;
			if (m<lLength)
				(*result) << ',';
			k = m-1;
		}
		else
		{
			conv = lData[k];
			(*result) << &conv;
			if (k<lLength-1)
				(*result) << ',';
		}
	}
	(*result).Finalize();
	return result;
}

//______________________________________________________________
// append & store operator
void _SimpleList::RequestSpace (long slots)
{
	if (slots>laLength)
	{
		laLength=(slots/MEMORYSTEP+1)*MEMORYSTEP;
		if (lData)
			checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
		else
			checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
	}
}
//______________________________________________________________
// append & store operator
void _List::InsertElement (BaseRef br, long insertAt, bool store)
{
	_SimpleList::InsertElement (br, insertAt, store);
}


//______________________________________________________________
// char* conversion

BaseRef _List::toStr(void)
{
	_String * s = new _String((unsigned long)20*(lLength+1),true);
	
	checkPointer (s);
	
	(*s)<<'{';
		
	for (unsigned long i = 0; (i<lLength); i++)
	{
		_String* t = (_String*)(((BaseRef*)lData)[i]->toStr());
		if (t)
		{
			(*s)<<t;
			DeleteObject (t);
		}
		if (i<lLength-1) 
			(*s)<<',';
	}
	(*s)<<'}';
	s->Finalize();
	return s;
}

//______________________________________________________________
// char* conversion

void _List::toFileStr(FILE* dest)
{
	fprintf (dest,"{");
		
	for (unsigned long i = 0; (i<lLength); i++)
	{
		((BaseRef*)lData)[i]->toFileStr(dest);
		if (i<lLength-1) fprintf (dest,",");
	}
	fprintf (dest,"}");
}

//______________________________________________________________
// char* conversion

BaseRef _SimpleList::toStr(void)
{
	if (lLength)
	{
		unsigned long ssi = _String::storageIncrement,
				 	  ma  = lLength*(1+log10((double)lLength));
		
		if ( ma > ssi)
			_String::storageIncrement = ma;
			
		_String * s = new _String (10L, true);
		
		checkPointer (s);
		
		(*s) << "{";
			
		for (unsigned long i = 0; i<lLength; i++)
		{
			char c[32];
			sprintf(c,"%ld",((long*)lData)[i]),
			(*s) << c;
			if (i<lLength-1) 
				(*s) << ',';
		}
		
		(*s) << '}';
		
		s->Finalize();
		_String::storageIncrement = ssi;		
		return s;
	}
	else
	{
		return new _String ("{}");
	}
}
	
//______________________________________________________________
 				
BaseRef _List::makeDynamic(void)
{
	_List * Res = new _List;
	checkPointer(Res);
	//lData = nil;
	memcpy ((char*)Res, (char*)this, sizeof (_List));
	Res->nInstances = 1;
	Res->lData = nil;
	Res->Duplicate (this);
	return Res;
}

//______________________________________________________________
 				
BaseRef _SimpleList::makeDynamic(void)
{
	_SimpleList * Res = new _SimpleList;
	checkPointer(Res);
	memcpy ((char*)Res, (char*)this, sizeof (_SimpleList));
	Res->nInstances = 1;
	Res->lData = nil;
	Res->Duplicate (this);
	return Res;
}
//______________________________________________________________

void	_List::bumpNInst (void)
{
	for (unsigned long i = 0; i<lLength; i++)
	{
		((BaseRef*)lData)[i]->nInstances++;
	}
}
 
//______________________________________________________________

long  _List::Find (BaseRef s, long startat)
{
	_String * st = (_String*)s;
	for (unsigned long i = startat; i<lLength; i++)
	{
		 _String * sp = (_String*)(((BaseRef*)lData)[i]->toStr());
		 
		 if (st->Equal(sp)) 
		 {
		 	DeleteObject(sp);
		 	return i;
		 }
		 DeleteObject(sp);
	}
	return -1;
}	

//______________________________________________________________

long  _List::FindString (BaseRef s, long startat, bool caseSensitive, long upTo)
{
	char * s1, *s2;
	long t = ((_String*)s)->sLength;
	
	if (upTo < 0 || upTo >= lLength)
		upTo = lLength-1;
	
	for (long i = startat; i<= upTo; i++)
	{
		 s1 = ((_String*)s)->sData;
		 if (((_String*)(((BaseRef*)lData)[i]))->sLength==t)
		 {
			 s2 = ((_String*)(((BaseRef*)lData)[i]))->sData;
			 long j = 0;
			 if (caseSensitive)
				 for (j=0; (*s1==*s2)&&(j<t); j++,s1++,s2++) ;
			 else
				 for (j=0; (toupper(*s1)==toupper(*s2))&&(j<t); j++,s1++,s2++) ;
			 
			 if (j==t)
			 	return i;
		  }
	}
	return -1;
}	
//______________________________________________________________

BaseRef  _List::Join (BaseRef spacer)
{
	_String *joined = new _String (256L,true);
	
	for (long k = 0; k < lLength; k++)
	{
		if (k)
			(*joined) << *(_String*)spacer;
		joined->AppendNewInstance((_String*) ((BaseRef*)lData)[k]->toStr());
	}
	
	joined->Finalize();
	return joined;
}

//______________________________________________________________

long  _List::BinaryFind (BaseRef s)
{
	_String * st = (_String*)s;
	long top=lLength-1, bottom=0, middle;
	
	if (top==-1) return -1;
	while (top>bottom)
	{
		middle = (top+bottom)/2;
		_String* stp = (_String*)(((BaseRef*)lData)[middle]->toStr());
		
		int		 cres = st->Compare (stp);
		DeleteObject (stp);
		
		if (cres < 0)
			top = middle==top?top-1:middle;
		else
			if (cres == 0)
				return middle;
			else
				bottom = middle==bottom?bottom+1:middle;
			
	}
	middle = top;
	_String* stp=(_String*)(((BaseRef*)lData)[middle]->toStr());
	if (st->Equal(stp))
	{
		DeleteObject(stp);
		return middle;
	}
	DeleteObject(stp);
	return -middle-2;
}	

//______________________________________________________________

long  _List::BinaryInsert (BaseRef s)
{
	if (!lLength)
	{
		InsertElement (s,0,true);
		return 0;
	}
		
	long pos = -BinaryFind (s)-2;
	if (pos<0) return -pos+2;
	_String *s1 = (_String*)s->toStr(), *s2 =(_String*) ((*this)(pos))->toStr();
	if (*s2<*s1)
	{
		pos++;
	}
	DeleteObject(s1);
	DeleteObject(s2);
	InsertElement (s,pos,true);
	return pos>=lLength?lLength-1:pos;
}

//______________________________________________________________

long		  _SimpleList::Min			(void)
{
	long res = LONG_MAX;
	for  (long e = 0; e < lLength; e++)
		if (lData[e] < res)
			res = lData[e];
	return res;
}

//______________________________________________________________

long		  _SimpleList::Max			(void)
{
	long res = LONG_MIN;
	for  (long e = 0; e < lLength; e++)
		if (lData[e] > res)
			res = lData[e];
	return res;
}

//______________________________________________________________

void		  _SimpleList::DebugVarList			(void)
{
	printf ("\nVariable list dump:\n");
	for  (long e = 0; e < lLength; e++)
	{
		if (lData[e] >= 0)
		{
			_Variable * theV = LocateVar (lData[e]);
			if (theV)
			{
				printf ("[%s]\n", theV->GetName()->getStr());
				continue;
			}
		}		
		printf ("[Empty]\n");
	}
}


//______________________________________________________________

_SimpleList*  _SimpleList::CountingSort (long upperBound, _SimpleList* ordering)
{
	if (ordering)
		ordering->Clear();
	
	if (lLength)
	{
		if (upperBound < 0)
			upperBound = Max()+1;

		_SimpleList buffer      (upperBound, 0, 0),
				  * result    =  new _SimpleList (lLength);
		for (long pass1 = 0; pass1 < lLength; pass1 ++)
			buffer.lData[lData[pass1]] ++;
		for (long pass2 = 1; pass2 < upperBound; pass2 ++)
			buffer.lData[pass2] += buffer.lData[pass2-1];
		if (ordering)
		{
			ordering->Populate (lLength, 0, 0);
			for (long pass3 = lLength-1; pass3 >=0; pass3--)
			{
				result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
				ordering->lData[buffer.lData[lData[pass3]]] = pass3;
			}
		}	
		else
			for (long pass3 = lLength-1; pass3 >=0; pass3--)
				result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
		result->lLength = lLength;

		return result;
	}
	return new _SimpleList;
}


//______________________________________________________________

long  _SimpleList::BinaryInsert (long n)
{
	if (!lLength)
	{
		(*this) << n;
		return 0;
	}
		
	long pos = -BinaryFind (n)-2;
	
	if (pos<0) 
		return -pos+2;
		
	if (lData[pos]<n)
		pos++;

	InsertElement ((BaseRef)n,pos,false,false);
	
	return pos>=lLength?lLength-1:pos;
}

//______________________________________________________________

long  _SimpleList::Find (long s, long startAt)
{
	for (unsigned long i = startAt; i<lLength; i++)
	{
		 if ( ((long*)(lData))[i] == s ) return i;
	}
	return -1;
}	

//______________________________________________________________

long  _SimpleList::FindStepping (long s, long step, long startAt)
{
	for (unsigned long i = startAt; i<lLength; i+=step)
		 if (lData[i] == s) 
		 	return i;
		 	
	return -1;
}	

//______________________________________________________________

void  _SimpleList::FilterRange (long lb, long ub)
{
	if (ub <= lb)
		Clear();
	else
	{
		_SimpleList toDelete;
		for (long k = 0; k < lLength; k++)
			if (lData[k] <= lb || lData[k] >= ub)
				toDelete << k;
		DeleteList (toDelete);
	}
}	

//______________________________________________________________

long  _SimpleList::BinaryFind (long s, long startAt)
{
	long top	=	lLength-1, 
		 bottom	=	startAt, 
		 middle;
	
	if (top == -1) 
		return -2;
	
	while (top>bottom)
	{
		middle = (top+bottom)/2;
		if (s<((long*)lData)[middle])
			top = middle==top?top-1:middle;
		else
			if (s>((long*)lData)[middle])
				bottom = middle==bottom?bottom+1:middle;
			else
				return middle;		
	}
	
	middle	   = top;
	long comp  = ((long*)lData)[middle]-s;
	if (!comp)
		return middle;
	
	return comp<0?-middle-3:-middle-2;
}	

//______________________________________________________________

void  _SimpleList::Sort (bool ascending)
{
	if (lLength<10) // use bubble sort
		BubbleSort();
	else
		QuickSort(0,lLength-1);
		
	if (!ascending)
	{
		long swap,i,j;
		for (i=0, j=lLength-1;i<j;i++,j--)
		{
			swap = ((long*)lData)[i];
			((long*)lData)[i]=((long*)lData)[j];
			((long*)lData)[j]=swap;
		}
	}
}	

	

//______________________________________________________________

void  _SimpleList::BubbleSort (void)
{
	bool done = false;
	long swap,i,j;
	while (!done)
	{
		done = true;
		for (i=lLength-1,j=i-1;i>0;i--,j--)
		{
			if (Compare(i,j)<0)
			{
				done = false;
				swap = ((long*)lData)[i];
				((long*)lData)[i]=((long*)lData)[j];
				((long*)lData)[j]=swap;
			}
		}
	}
}	

//______________________________________________________________

void  _SimpleList::QuickSort (long from, long to)
{
	long middle = (from+to)/2, 
		 middleV = ((long*)lData)[middle], 
		 top = to,
		 bottommove = 1, 
		 topmove = 1, 
		 temp,
		 i;
		 
	if (middle)
		//while ((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
		while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>0))
		{
			bottommove++;
		}
		
	if (from<to)
		//while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
		while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<0))
		{
			topmove++;
		}
	// now shuffle
	for (i=from; i<middle-bottommove; i++)
	{
		if (Compare(i,middle)>0)
		{
			temp = ((long*)lData)[middle-bottommove];
			((long*)lData)[middle-bottommove] = ((long*)lData)[i];
			((long*)lData)[i]=temp;
			bottommove++;
			
			//while ((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
			while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>0))
			{
				bottommove++;
			}
		}
	}
	
	for (i=middle+topmove+1; i<=top; i++)
	{
		if (Compare(i,middle)<0)
		{
			temp = ((long*)lData)[middle+topmove];
			((long*)lData)[middle+topmove] = ((long*)lData)[i];
			((long*)lData)[i]=temp;
			topmove++;
			
			//while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
			while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<0))
			{
				topmove++;
			}
		}
	}
	
	if (topmove==bottommove)
	{
		for (i=1; i<bottommove; i++)
		{
			temp = ((long*)lData)[middle+i];
			((long*)lData)[middle+i] = ((long*)lData)[middle-i];
			((long*)lData)[middle-i]=temp;
		}
	}
	else
	if (topmove>bottommove)
	{
		long shift = topmove-bottommove;
		for (i=1; i<bottommove; i++)
		{
			temp = ((long*)lData)[middle+i+shift];
			((long*)lData)[middle+i+shift] = ((long*)lData)[middle-i];
			((long*)lData)[middle-i]=temp;
		}
		for (i=0; i<shift; i++)
		{
			((long*)lData)[middle+i]=((long*)lData)[middle+i+1];
		}
		middle+=shift;
		((long*)lData)[middle]=middleV;
	}
	else
	{
		long shift = bottommove-topmove;
		for (i=1; i<topmove; i++)
		{
			temp = ((long*)lData)[middle-i-shift];
			((long*)lData)[middle-i-shift] = ((long*)lData)[middle+i];
			((long*)lData)[middle+i]=temp;
		}
		for (i=0; i<shift; i++)
		{
			((long*)lData)[middle-i]=((long*)lData)[middle-i-1];
		}
		middle-=shift;
		((long*)lData)[middle]=middleV;
	}
	if (to>middle+1)
		QuickSort (middle+1,top);
	if (from<middle-1)
		QuickSort (from,middle-1);
}
//______________________________________________________________

long  _SimpleList::Compare (long i, long j)
{
	long    v1 = ((long*)lData)[i],
			v2 = ((long*)lData)[j];
			

	if (v1<v2)
		return -1;
	else
		if (v1==v2) 
			return 0;
		else
			return 1;
	
		
	//return ((long*)lData)[i]-((long*)lData)[j];
}

//______________________________________________________________

long  _SimpleList::Compare (BaseRef i, long j)
{
	long    v1 = (long)i,
			v2 = ((long*)lData)[j];
			

	if (v1<v2)
		return -1;
	else
		if (v1==v2) 
			return 0;
		else
			return 1;
	
	//return (long)i-((long*)lData)[j];
}

//______________________________________________________________

long  _List::Compare (long i, long j)
{
	_String				*si = (_String*)lData[i], 
						*sj = (_String*)lData[j];
						
	return	si->Compare(sj);
}

//______________________________________________________________

long  _List::Compare (BaseRef i, long j)
{
	_String				*sj = (_String*)lData[j], 
						*si = (_String*)i;
						
	return	si->Compare(sj);
}

//______________________________________________________________

long  _List::FreeUpMemory (long requestedBytes)
{
	long freed = 0;
	for (unsigned long i = 0; i<lLength; i++)
	{
		BaseRef t = ((BaseRef*)lData)[i];
		freed+=t->FreeUpMemory(requestedBytes-freed);
		if (freed>=requestedBytes)
		{
			return freed;
		}
	}
	return freed;
}


//______________________________________________________________

void  _List::Clear (bool completeClear)
{
	if (nInstances<=1)
	{
		for (unsigned long i = 0; i<lLength; i++)
			DeleteObject (((BaseRef*)lData)[i]);
		_SimpleList::Clear(completeClear);
			
	}
	else
		nInstances--;	
}

//______________________________________________________________

void  _SimpleList::Clear (bool completeClear)
{
	if (nInstances<=1)
	{
		lLength = 0;
		if (completeClear)
		{
			laLength = 0;
			if (lData)
				free (lData);
			lData = nil;
		}
	}
	else
		nInstances--;	
}

//______________________________________________________________

void  _List::Delete (long index)
//delete item at index (>=0)
{
	if ((index>=0)&&(index<lLength))
	{
		BaseRef theObj = ((BaseRef*)lData)[index];
		DeleteObject (theObj);
		lLength--;
		if (lLength-index)
			for (unsigned long i = index; i < lLength; i++)
				lData[i] = lData[i+1];
			//memcpy ((Ptr)lData+sizeof(BaseRef)*(index),(Ptr)lData+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
	}
	if (laLength-lLength>MEMORYSTEP)
	{	
		laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
		lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
	}
		
}

//______________________________________________________________

void  _List::DeleteList (const _SimpleList& toDelete)
//delete item at index (>=0)
{
	if (toDelete.lLength)
	{
		long k = 0;
		for (long i = 0; i<lLength; i++)
		{
			if (k<toDelete.lLength && i==toDelete.lData[k])
			{
				DeleteObject (((BaseRef*)lData)[i]);
				//if (k<toDelete.lLength)
				k++;
			}
			else
				((BaseRef*)lData)[i-k] = ((BaseRef*)lData)[i];
		}
		lLength -= toDelete.lLength;
		if (laLength-lLength>MEMORYSTEP)
		{	
			laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
			lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
		}	
	}	
}


//______________________________________________________________

void  _List::Replace (long index, BaseRef newObj, bool dup)
{
	if ((index>=0)&&(index<lLength))
	{
		BaseRef theObj = ((BaseRef*)lData)[index];
		DeleteObject (theObj);
		((BaseRef*)lData)[index] = dup?newObj->makeDynamic():newObj;
	}	
}

//______________________________________________________________

void  _SimpleList::Delete (long index, bool compact)
//delete item at index (>=0)
{
	if (index>=0 && index<lLength)
	{
		lLength--;
		if (lLength-index)
			memmove ((Ptr)lData+sizeof(BaseRef)*(index),(Ptr)lData+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
	}
	if (compact && laLength-lLength>MEMORYSTEP)
	{	
		laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
		if (laLength)
			lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
		else
		{
			free (lData);
			lData = nil;
		}
	}
		
}

//______________________________________________________________

void  _SimpleList::TrimMemory (void)
//delete item at index (>=0)
{
	if (laLength>lLength)
	{	
		laLength = lLength;
		if (laLength)
		{
			if (lData)
				lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
			else
				lData = (long*)MemAllocate (laLength*sizeof(Ptr));
			if (!lData)
				checkPointer (lData);
		}		
		else
		{
			if (lData)
			{
				free (lData);
				lData = nil;
			}
		}
	}	
}

//______________________________________________________________

void  _SimpleList::DeleteDuplicates (void)
// delete duplicates from a sorted list
{
	if (lLength>1)
	{
		_SimpleList noDups;
		
		long 	lastValue = lData[0]+1;
		for (long k=0; k<lLength; k++)
		{
			long thisValue = lData[k];
			if (thisValue!=lastValue)
			{
				noDups << thisValue;
				lastValue = thisValue;
			}
		}
		
		if (noDups.lLength!=lLength)
			Duplicate (&noDups);
	}
		
}

//______________________________________________________________

void  _SimpleList::DeleteList (const _SimpleList& toDelete)
//delete items from a sorted list
{
	if (toDelete.lLength)
	{
		long k = 0;
		for (long i = 0; i<lLength; i++)
		{
			if (k<toDelete.lLength && i==toDelete.lData[k])
				//if (k<toDelete.lLength)
				k++;
			else
				lData[i-k] = lData[i];
		}
		lLength -= toDelete.lLength;
	}	
	if (laLength-lLength>MEMORYSTEP)
	{	
		laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
		if (laLength)
			lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
		else
		{
			free (lData);
			lData = nil;
		}
	}	
}
//______________________________________________________________

void  _SimpleList::Displace (long start, long end, long delta)
//shift the range from start to end
{
	if (start<0) 
		start = 0;
	else
		if (start>=lLength) 
			start = lLength-1;
			
	if (end<0) 
		end = lLength-1;
	else
		if (end>=lLength) 
			end = lLength-1;
	
	if ((end-start>=0)&&delta&&(end-start<lLength-1))
	// stuff to do
	{
		if (delta>0) // shift up
		{
			if (lLength-end<=delta)
				delta = lLength-end-1;
		}
		else
		{
			if (start-delta<0)
			{
				delta = start;
			}
		}
		if (delta)
		{
			long i,j,delta2 = end-start+1;
			_SimpleList swapList ((unsigned long)(end-start+1));
			for (i=start; i<=end; i++)
				swapList << lData[i];
			if (delta>0)
				for (i=end+1; i<=end+delta; i++)
					lData[i-delta2] = lData[i];
			else
				for (i=start-1; i>=start+delta; i--)
					lData[i+delta2] = lData[i];
			for (i=start+delta, j=0;i<=end+delta; i++,j++)
				lData[i] = swapList.lData[j];
		}		
		
	}		
}
//______________________________________________________________

void  _SimpleList::PermuteWithReplacement (long blockLength)
// create a permutation of the list's elements with possible repetitions
{
	unsigned long blockCount = lLength/blockLength;
	_SimpleList	  result ((unsigned long)(blockCount*blockLength));
	if (blockLength>1)
		for (long i = 0; i<blockCount; i++)
		{
			unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
			sample *= blockLength;
			for (long j = 0; j<blockLength; j++,sample++)
				result<<lData[sample];
		}
	else
		for (long i = 0; i<blockCount; i++)
		{
			unsigned long sample = genrand_real2()*blockCount;
			result<<lData[sample];
		}

	Clear();
	Duplicate(&result);
		
}

//______________________________________________________________

void  _SimpleList::Permute (long blockLength)
// create a permutation of the list's elements
{
	unsigned long blockCount = lLength/blockLength;
	
	if (blockLength>1)
	{
		/*_SimpleList	  result ((unsigned long)(blockCount*blockLength));
		while (blockCount)
		{
			unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
			sample *= blockLength;
			for (long j = 0; j<blockLength; j++)
			{
				result<<lData[sample];
				Delete(sample);
			}
			blockCount --;
		}
		Duplicate(&result);	*/	

		for (unsigned long k=0; k<blockCount-1; k=k+1)
		{
			unsigned long k2 = genrand_real2()*(blockCount-k);
			if (k2)
			{
				k2 += k;
				k2 *= blockLength;

				for (long j = 0; j<blockLength; j++)
				{
					long t = lData[k2+j];
					lData[k2+j] = lData[k*blockLength+j];
					lData[k*blockLength+j] = t;
				}
			}
		}

	}
	else
	{		
		for (unsigned long k=0; k<blockCount-1; k=k+1)
		{
			unsigned long k2 = genrand_real2()*(blockCount-k);
			if (k2)
			{
				k2+=k;
				long t = lData[k2];
				lData[k2] = lData[k];
				lData[k] = t;
			}
		}
		/*{
			while (blockCount)
			{
				unsigned long sample = genrand_real2()*blockCount;
				result<<lData[sample];
				Delete(sample);
				blockCount --;
			}
		}*/
	}
}


//______________________________________________________________

bool  _SimpleList::NChooseKInit (_SimpleList& state, _SimpleList& store, unsigned long stride, bool algorithm)
// together with the next function
// implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
{
	if (stride <= lLength && lLength)
	{
		state.Clear();
		state.RequestSpace (stride+3);
		state << stride;
		store.Clear();
		store.RequestSpace (stride);
		return true;
	}
	return false;
}

//______________________________________________________________

bool  _SimpleList::NChooseK (_SimpleList& state, _SimpleList& store)
// implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
{
	if (state.lLength == 1) // first pass
	{
		state << 0;              // m
		state << state.lData[0]; // h
		state.lLength = state.lData[0]+3;
		store.lLength = state.lData[0];
		if (store.lLength == 0)
			return false;
	}
	else
	{
		if (state.lData[1] < lLength - state.lData[2])
			state.lData[2] = 0;
		state.lData[2] ++; 
		state.lData[1] = state.lData[3 + state.lData[0] - state.lData[2]] + 1;
	}	
	for (long j=1; j <= state.lData[2]; j++)
	{
		long  anIndex	= j+state.lData[0]-state.lData[2],
			  anIndex2	= state.lData[1]+j;
		state.lData[anIndex+2]   = anIndex2-1;
		store.lData[anIndex-1]   = lData[anIndex2-1];
	}	
	return state.lData[3] < lLength-state.lData[0];
}

//______________________________________________________________
void	_SimpleList::Swap (long i, long j)
{
	if ( i>=lLength || j>=lLength )
		return;	

	void * pt = ((void**)lData)[j];
	((void**)lData)[j] = ((void**)lData)[i];
	((void**)lData)[i] = pt;
}

//______________________________________________________________
void	_SimpleList::Flip ()
{
	for (long k=0, l=lLength-1; k<l; k++,l--)
	{
		void * pt = ((void**)lData)[k];
		((void**)lData)[k] = ((void**)lData)[l];
		((void**)lData)[l] = pt;
	}
}
//______________________________________________________________
 				 			
void	    SortLists (_SimpleList* ref, _SimpleList* index)
{
	if ((*ref).lLength!=index->lLength) return;
	if ((*ref).lLength<=10)
	{
		bool done = false;
		
		while (!done)
		{
			done = true;
			for (long i=1; i<(*ref).lLength; i++)
			{
				if (ref->Compare(i-1,i)>0)
				{
					long swap;
					swap = ((long*)ref->lData)[i];
					((long*)ref->lData)[i]=((long*)ref->lData)[i-1];
					((long*)ref->lData)[i-1]=swap;
					swap = ((long*)index->lData)[i];
					((long*)index->lData)[i]=((long*)index->lData)[i-1];
					((long*)index->lData)[i-1]=swap;
					done = false;
				}
			}
		}
	}
	else
	{
		(*ref).RecursiveIndexSort (0, (*ref).lLength-1,index);
	}
}

//______________________________________________________________

void	_SimpleList::RecursiveIndexSort (long from, long to, _SimpleList* index)
{
	long middle = (from+to)/2, middleV = lData[middle],
		 bottommove = 1, topmove = 1, temp,i, imiddleV = (*index)(middle);
	long *idata = (*index).quickArrayAccess();
	if (middle)
		while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>=0))
			bottommove++;

	if (from<to)
		while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<=0))
			topmove++;

	// now shuffle
	for (i=from; i<middle-bottommove; i++)
	{
		if (Compare(i,middle)>=0)
		{
			temp = lData[middle-bottommove];
			lData[middle-bottommove] = lData[i];
			lData[i]=temp;
			temp = idata[middle-bottommove];
			idata[middle-bottommove] = idata[i];
			idata[i]=temp;
			bottommove++;
			while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>=0))
			{
				bottommove++;
			}
		}
	}
	
	for (i=middle+topmove+1; i<=to; i++)
	{
		if (Compare(i,middle)<=0)
		{
			temp = lData[middle+topmove];
			lData[middle+topmove] = lData[i];
			lData[i]=temp;
			temp = idata[middle+topmove];
			idata[middle+topmove] = idata[i];
			idata[i]=temp;
			topmove++;
			while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<=0))
			{
				topmove++;
			}
		}
	}
	
	if (topmove==bottommove)
	{
		for (i=1; i<bottommove; i++)
		{
			temp = lData[middle+i];
			lData[middle+i] = lData[middle-i];
			lData[middle-i]=temp;
			temp = idata[middle+i];
			idata[middle+i] = idata[middle-i];
			idata[middle-i]=temp;
		}
	}
	else
	if (topmove>bottommove)
	{
		long shift = topmove-bottommove;
		for (i=1; i<bottommove; i++)
		{
			temp = lData[middle+i+shift];
			lData[middle+i+shift] = lData[middle-i];
			lData[middle-i]=temp;
			temp = idata[middle+i+shift];
			idata[middle+i+shift] = idata[middle-i];
			idata[middle-i]=temp;
		}
		for (i=0; i<shift; i++)
		{
			lData[middle+i]=lData[middle+i+1];
			idata[middle+i]=idata[middle+i+1];
		}
		middle+=shift;
		lData[middle]=middleV;
		idata[middle]=imiddleV;
	}
	else
	{
		long shift = bottommove-topmove;
		for (i=1; i<topmove; i++)
		{
			temp = lData[middle-i-shift];
			lData[middle-i-shift] = lData[middle+i];
			lData[middle+i]=temp;
			temp = idata[middle-i-shift];
			idata[middle-i-shift] = idata[middle+i];
			idata[middle+i]=temp;
		}
		for (i=0; i<shift; i++)
		{
			lData[middle-i]=lData[middle-i-1];
			idata[middle-i]=idata[middle-i-1];
		}
		middle-=shift;
		lData[middle]=middleV;
		idata[middle]=imiddleV;
	}
	if (to>middle+1)
		RecursiveIndexSort (middle+1,to, index);
	if (from<middle-1)
		RecursiveIndexSort (from,middle-1, index);
}

//______________________________________________________________

void	_SimpleList::Union (_SimpleList& l1, _SimpleList& l2)
// compute the union of two sorted lists
// each repeat appears exactly once
{
	if (lLength)
		Clear();
		
	long  c1 = 0,
		  c2 = 0;
		  
	while (c1<l1.lLength && c2<l2.lLength)
	{
		while (l1.lData[c1]<l2.lData[c2])
		{
			(*this) << l1.lData[c1++];
			if (c1==l1.lLength)
				break;
		}
		
		if (c1==l1.lLength)
			break;
		
		while (l1.lData[c1]==l2.lData[c2]) 
		{
			(*this) << l1.lData[c1++];
			c2++;
			if (c1==l1.lLength || c2==l2.lLength)
				break;
		}
		
		if (c1==l1.lLength || c2==l2.lLength)
			break;
		
		while (l2.lData[c2]<l1.lData[c1])
		{
			(*this) << l2.lData[c2++];
			if (c2==l2.lLength)
				break;
		}
	}
	
	while (c1<l1.lLength)
		(*this) << l1.lData[c1++];
	while (c2<l2.lLength)
		(*this) << l2.lData[c2++];
	
}

//______________________________________________________________

void	_SimpleList::Intersect (_SimpleList& l1, _SimpleList& l2)
// compute the intersection of two sorted lists
{
	if (lLength)
		Clear();
		
	long  c1 = 0,
		  c2 = 0;
		  
	while (c1<l1.lLength && c2<l2.lLength)
	{
		while (l1.lData[c1]<l2.lData[c2]) 
		{
			c1++;
			if (c1==l1.lLength)
				break;
		}
		if (c1==l1.lLength)
			break;
		
		while (l1.lData[c1]==l2.lData[c2]) 
		{
			(*this) << l1.lData[c1++];
			c2++;
			if (c1==l1.lLength || c2==l2.lLength)
				break;
		}
		if (c1==l1.lLength || c2==l2.lLength)
			break;
		while (l2.lData[c2]<l1.lData[c1])
		{
			c2++;
			if (c2==l2.lLength)
				break;
		}
	}
}

//______________________________________________________________

long	_SimpleList::CountCommonElements (_SimpleList& l1, bool yesNo)
// compute the number of shared of two sorted lists
{
	long  c1 	= 0,
		  c2 	= 0,
		  res 	= 0;
		  
	
	while (c1<l1.lLength && c2<lLength)
	{
		while (l1.lData[c1]<lData[c2]) 
		{
			c1++;
			if (c1==l1.lLength)
				break;
		}
		if (c1==l1.lLength)
			break;
		
		while (l1.lData[c1]==lData[c2]) 
		{
			c2++;
			if (yesNo)
				return 1;
			else
				res++;
			if (c1==l1.lLength || c2==lLength)
				break;
		}
		if (c1==l1.lLength || c2==lLength)
			break;
		while (lData[c2]<l1.lData[c1])
		{
			c2++;
			if (c2==lLength)
				break;
		}
	}

	return res;
}

//______________________________________________________________

void	_List::Intersect (_List& l1, _List& l2, _SimpleList* idx, _SimpleList* idx2)
// compute the union of two sorted lists
// each repeat appears exactly once
{
	if (lLength)
		Clear();
		
	long  c1 = 0,
		  c2 = 0;
		  
	while (c1<l1.lLength && c2<l2.lLength)
	{
		while (c1<l1.lLength && ((_String*)l1(c1))->Compare((_String*)l2(c2))<0) 
			c1++;
		if (c1==l1.lLength)
			break;
		
		while (c1<l1.lLength && c2<l2.lLength &&((_String*)l1(c1))->Equal((_String*)l2(c2))) 
		{
			if (idx)
				(*idx) << c1;
			if (idx2)
				(*idx2) << c2;
				
			(*this) << l1(c1++);
			c2++;
		}
		if (c1==l1.lLength || c2==l2.lLength)
			break;
		while (c2<l2.lLength && ((_String*)l2(c2))->Compare((_String*)l1(c1))<0)
			c2++;
	}
}

//______________________________________________________________

void	_SimpleList::XOR (_SimpleList& l1, _SimpleList& l2)
// compute the union of two sorted lists
// each repeat appears exactly once
{
	if (lLength)
		Clear();
		
	long  c1 = 0,
		  c2 = 0;
		  
	while ((c1<l1.lLength)&&(c2<l2.lLength))
	{
		while (c1<l1.lLength && l1.lData[c1]<l2.lData[c2])
			(*this) << l1.lData[c1++];
		if (c1==l1.lLength)
			break;
		while (c1<l1.lLength && c2<l2.lLength && l1.lData[c1]==l2.lData[c2] ) 
		{
			c1++;
			c2++;
		}
		if (c1==l1.lLength||c2==l2.lLength)
			break;
		
		while (c2<l2.lLength && l2.lData[c2]<l1.lData[c1]) 
			(*this) << l2.lData[c2++];
	}
	
	while (c1<l1.lLength)
		(*this) << l1.lData[c1++];
	while (c2<l2.lLength)
		(*this) << l2.lData[c2++];
}

//______________________________________________________________

void	_SimpleList::Subtract (_SimpleList& l1, _SimpleList& l2)
// compute the union of two sorted lists
// each repeat appears exactly once
{
	if (lLength)
		Clear();
		
	long  c1 = 0,
		  c2 = 0;
		  
	while (c1<l1.lLength && c2<l2.lLength)
	{
		while (c1<l1.lLength && l1.lData[c1]<l2.lData[c2]) 
			(*this) << l1.lData[c1++];
		if (c1==l1.lLength)
			break;
		while ( c1<l1.lLength && c2<l2.lLength  && l1.lData[c1]==l2.lData[c2] ) 
		{
			c1++;
			c2++;
		}
		if (c1==l1.lLength || c2==l2.lLength)
			break;
		while (c2<l2.lLength && l2.lData[c2]<l1.lData[c1]) 
			c2++;
	}
	
	while (c1<l1.lLength)
		(*this) << l1.lData[c1++];
}


//______________________________________________________________
// AVL Lists
//______________________________________________________________

_AVLList::_AVLList (_SimpleList* d)
{
	dataList = d;
	root	 = -1;
}

//______________________________________________________________

_AVLListX::_AVLListX (_SimpleList* d):_AVLList(d)
{
}

//______________________________________________________________

_AVLListXL::_AVLListXL (_SimpleList* d):_AVLList(d)
{
}

//______________________________________________________________

long	_AVLListX::GetXtra (long d)
{
	return xtraD.lData[d];
}

//______________________________________________________________

BaseRef	_AVLListXL::GetXtra (long d)
{
	return xtraD(d);
}

//______________________________________________________________

void	_AVLListX::SetXtra (long i, long d)
{
	xtraD.lData[i] = d;
}

//______________________________________________________________

void	_AVLListXL::SetXtra (long i, BaseRef d, bool dup)
{
	xtraD.Replace (i,d, dup);
}

//______________________________________________________________

long  _AVLList::Find (BaseRef obj)
{
	long curNode = root;
		 
	while (curNode>=0)
	{
		long comp = dataList->Compare (obj,curNode);
		
		if (comp<0)
			curNode = leftChild.lData[curNode];
		else
			if (comp>0)
				curNode = rightChild.lData[curNode];
			else
				return curNode;
	}

	return -1;
}

//______________________________________________________________

long  _AVLList::FindLong (long obj)
{
	long curNode = root;
	
	while (curNode>=0)
	{
		long comp = dataList->lData[curNode];
		
		if (obj<comp)
			curNode = leftChild.lData[curNode];
		else
			if (obj>comp)
				curNode = rightChild.lData[curNode];
			else
				return curNode;
	}
	
	return -1;
}

//______________________________________________________________

char  _AVLList::FindBest (BaseRef obj, long& lastNode)
{
	long curNode  = root,
		 comp     = 1;
			 
	while (curNode>=0 && comp)
	{
		comp = dataList->Compare (obj,curNode);
		lastNode = curNode;
		
		if (comp<0)
			curNode = leftChild.lData[curNode];
		else
			if (comp>0)
				curNode = rightChild.lData[curNode];
			else
				return 0;
	}
	
	return comp;
}

//______________________________________________________________

long  _AVLList::Find (BaseRef obj, _SimpleList& hist)
{
	long curNode = root;
		 
	while (curNode>=0)
	{
		long comp = dataList->Compare (obj,curNode);
		
		if (comp<0)
		{
			hist << curNode;
			curNode = leftChild.lData[curNode];
		}
		else
			if (comp>0)
			{
				hist << curNode;
				curNode = rightChild.lData[curNode];
			}
			else
				return curNode;
	}

	return -1;
}

//______________________________________________________________

long  _AVLList::Next (long d, _SimpleList& hist)
{
	if (d >= 0)
	{
		if (rightChild.lData [d] >= 0)
		{
			hist << d;
			d = rightChild.lData [d];
			while (leftChild.lData[d] >= 0)
			{
				hist << d;
				d = leftChild.lData[d]; 
			}
			return d;
		}
		else
		{
			while (hist.countitems())
			{
				long x = hist.lData[hist.lLength-1];
				
				
				hist.Delete (hist.lLength-1);
				
				if (rightChild.lData[x] != d)
					return x;
					
				d = x;
			}
			
			return -1;
		}
	}
	
	d = root;
	while (d >= 0 && leftChild.lData[d] >=0)
		d = leftChild.lData[d];
		
	return d;	
}

//______________________________________________________________

long  _AVLList::First (void)
{
	long   d = root;
	while (d >= 0 && leftChild.lData[d] >=0)
		d = leftChild.lData[d];
		
	return d;	
}

//______________________________________________________________

long  _AVLList::Last (void)
{
	long   d = root;
	while (d >= 0 && rightChild.lData[d] >=0)
		d = rightChild.lData[d];
		
	return d;	
}

//______________________________________________________________

long  _AVLList::Prev (long d, _SimpleList& hist)
{
	if (d >= 0)
	{
		if (leftChild.lData [d] >= 0)
		{
			hist << d;
			d = leftChild.lData [d];
			while (rightChild.lData[d] >= 0)
			{
				hist << d;
				d = rightChild.lData[d]; 
			}
			return d;
		}
		else
		{
			while (hist.countitems())
			{
				long x = hist.lData[hist.lLength-1];
				
				hist.Delete (hist.lLength-1);
				
				if (leftChild.lData[x] != d)
					return x;
					
				d = x;
			}
			
			return -1;
		}
	}
	
	d = root;
	while (d >= 0 && rightChild.lData[d] >=0)
		d = rightChild.lData[d];
		
	return d;
	
}

//______________________________________________________________

void  _AVLList::ReorderList (_SimpleList *s)
{
	_SimpleList reorderMe ((unsigned long)(dataList->lLength-emptySlots.lLength+1)),
				nodeStack ((unsigned long)32);
	
	long		curNode = root;
	
	while (1)
	{
		while (curNode >= 0)
		{
			nodeStack << curNode;
			curNode = leftChild.lData[curNode];
		}
		if (long h = nodeStack.lLength)
		{
			h--;
			curNode = nodeStack.lData[h];
			if (s)
				(*s) << curNode;
			reorderMe.InsertElement (((BaseRef*)dataList->lData)[curNode],-1,false,false);
			curNode = rightChild.lData[curNode];
			nodeStack.Delete (h, false);
		}	
		else
			break;
	}
	
	reorderMe.TrimMemory ();
	
	long* t 			= dataList->lData;
	dataList->lData 	= reorderMe.lData;
	dataList->lLength 	= reorderMe.lLength;
	dataList->laLength 	= reorderMe.laLength;
	reorderMe.lData 	= t;
}

//______________________________________________________________

void  _AVLList::ConsistencyCheck (void)
{
	_SimpleList nodeStack ((unsigned long)32);
	
	long		curNode  = root,
				lastNode = -1,
				checkCount = 0;
	
	while (1)
	{
		while (curNode >= 0)
		{
			nodeStack << curNode;
			curNode = leftChild.lData[curNode];
			if (curNode >= (long)dataList->lLength)
			{
				WarnError ("Failed Constistency Check in _AVLList");
				return;
			}
			
		}
		if (long h = nodeStack.lLength)
		{
			if (h>3*log (1.+countitems()))
			{
				WarnError ("Failed Constistency Check in _AVLList");
				return;
			}
			h--;
			curNode = nodeStack.lData[h];
			if (lastNode >= 0 && curNode >= 0)
			{
				if (dataList->Compare (Retrieve (lastNode), curNode) >= 0)
				{
					WarnError ("Failed Constistency Check in _AVLList");
					return;
				}
				checkCount++;
			}
			if ((balanceFactor.lData[curNode] < -1)||(balanceFactor.lData[curNode] > 1))
			{
				WarnError ("Failed Constistency Check in _AVLList");
				return;
			}
			lastNode = curNode;
			curNode = rightChild.lData[curNode];
			if (curNode >= (long)dataList->lLength)
			{
				WarnError ("Failed Constistency Check in _AVLList");
				return;
			}
			nodeStack.Delete (h, false);
		}	
		else
			break;
	}
	
	if (dataList->lLength && (dataList->lLength > checkCount + 1 + emptySlots.lLength))
	{
		WarnError ("Failed Constistency Check in _AVLList");
		return;
	}
		
}

//______________________________________________________________

long  _AVLList::Traverser (_SimpleList &nodeStack, long& t, long r)
{
	if 	(r >= 0)
	{
		 t = r;
		 nodeStack.Clear();
	}
	
	while (t >= 0)
	{
		nodeStack << t;
		t = leftChild.lData[t];
	}
	if (long h = nodeStack.lLength)
	{
		h--;
		t = nodeStack.lData[h];
		r = t;
		t = rightChild.lData[t];
		nodeStack.Delete (h, false);
		return r;
	}	
	
	return -1;
}

//______________________________________________________________

BaseRef  _AVLList::toStr (void)
{
	_String * str = new _String (128L, true);
	checkPointer (str);
	
	if (countitems() == 0)
		(*str) << "Empty Associative List";
	else
	{
		_SimpleList	 hist;
		long		 ls, cn;
		
		cn = Traverser (hist,ls,root);
		
		while (cn>=0)
		{
			long keyVal = (long)Retrieve (cn);
			(*str) << _String(keyVal);
			(*str) << '\n';
			cn = Traverser (hist,ls);
		}
	}
	
	str->Finalize();
	return str;
}

//______________________________________________________________

BaseRef _AVLList::Retrieve (long idx)
{
	return ((BaseRef*)dataList->lData)[idx];
}

//______________________________________________________________

void _AVLList::Clear (bool cL)
{
	if (cL)
		((_List*)dataList)->Clear();
	else
		dataList->Clear();
		
	emptySlots.Clear();
	root = -1;
	leftChild.Clear();
	rightChild.Clear();
	balanceFactor.Clear();
}

//______________________________________________________________

void _AVLListX::Clear (bool cL)
{
	xtraD.Clear();
	_AVLList::Clear(cL);
}

//______________________________________________________________

BaseRef _AVLListX::toStr (void)
{
	_String * str = new _String (128L, true);
	checkPointer (str);
	
	if (countitems() == 0)
		(*str) << "Empty Associative List";
	else
	{
		_SimpleList	 hist;
		long		 ls, cn;
		
		cn = Traverser (hist,ls,root);
		
		while (cn>=0)
		{
			_String * keyVal = (_String*)Retrieve (cn);
			(*str) << keyVal;
			(*str) << " : ";
			(*str) << _String(GetXtra (cn));
			(*str) << '\n';
			cn = Traverser (hist,ls);
		}
	}
	
	str->Finalize();
	return str;
}

//______________________________________________________________

BaseRef _AVLListXL::toStr (void)
{
	_String * str = new _String (128L, true);
	checkPointer (str);
	
	if (countitems() == 0)
		(*str) << "Empty Associative List";
	else
	{
		_SimpleList	 hist;
		long		 ls, cn;
		
		cn = Traverser (hist,ls,root);
		
		while (cn>=0)
		{
			_String * keyVal = (_String*)Retrieve (cn);
			(*str) << keyVal;
			(*str) << " : ";
			(*str) << (_String*)GetXtra (cn);
			(*str) << '\n';
			cn = Traverser (hist,ls);
		}
	}
	
	str->Finalize();
	return str;
}


//______________________________________________________________

void _AVLListXL::Clear (bool cL)
{
	xtraD.Clear();
	_AVLList::Clear(cL);
}

//______________________________________________________________

long  _AVLList::InsertData (BaseRef b, long, bool)
{
	long w = (long)emptySlots.lLength - 1,
		 n;
		 
	if (w>=0)
	{
		n = emptySlots.lData[w]; 
		emptySlots.Delete (w);
		leftChild.lData[n] = -1;
		rightChild.lData[n] = -1;
		balanceFactor.lData[n] = 0;
		((BaseRef*)dataList->lData)[n] = b;
	}
	else	
	{
		n = dataList->lLength;
		dataList->InsertElement (b,-1,false,false);
		leftChild  << -1;
		rightChild << -1;
		balanceFactor << 0;
	}
	return n;
}

//______________________________________________________________

long  _AVLListX::InsertData (BaseRef b, long d, bool)
{
	long w = (long)emptySlots.lLength - 1,
		 n;
		 
	if (w>=0)
	{
		n = emptySlots.lData[w]; 
		emptySlots.Delete (w);
		leftChild.lData[n] = -1;
		rightChild.lData[n] = -1;
		balanceFactor.lData[n] = 0;
		xtraD.lData[n] = d;
		((BaseRef*)dataList->lData)[n] = b;
	}
	else	
	{
		n = dataList->lLength;
		dataList->InsertElement (b,-1,false,false);
		leftChild  << -1;
		rightChild << -1;
		balanceFactor << 0;
		xtraD << d;
	}
	return n;
}

//______________________________________________________________

long  _AVLListXL::InsertData (BaseRef b, long xl, bool cp)
{
	long w = (long)emptySlots.lLength - 1,
		 n;
		 
	BaseRef x = (BaseRef)xl;
		 
	if (w>=0)
	{
		n = emptySlots.lData[w]; 
		emptySlots.Delete (w);
		leftChild.lData[n] = -1;
		rightChild.lData[n] = -1;
		balanceFactor.lData[n] = 0;
		((BaseRef*)xtraD.lData)[n] = x;
		if (cp)
			x->nInstances++;
		((BaseRef*)dataList->lData)[n] = b;
	}
	else	
	{
		n = dataList->lLength;
		dataList->InsertElement (b,-1,false,false);
		leftChild  << -1;
		rightChild << -1;
		balanceFactor << 0;
		xtraD << x;
		if (!cp)
			x->nInstances--;
	}
	return n;
}

//______________________________________________________________

unsigned long _AVLList::countitems (void)
{
	return dataList->lLength - emptySlots.lLength;
}

//______________________________________________________________

void _AVLListXL::DeleteXtra (long i)
{
	DeleteObject (((BaseRef*)xtraD.lData)[i]);
	(((BaseRef*)xtraD.lData)[i]) = nil;
}

//______________________________________________________________

void _AVLListX::DeleteXtra (long i)
{
	xtraD.lData[i] = -1;
}

//______________________________________________________________

void _AVLListX::PopulateFromList (_List& src)
{
	Clear(true);
	for (long k = 0; k < src.lLength; k++)
		Insert (src(k)->makeDynamic(),k,false);
}


//______________________________________________________________

long  _AVLList::Insert (BaseRef b, long xtra,bool cp,bool clear)
{
	if (dataList->lLength-emptySlots.lLength)
	{
		long		y = root,
					z = -1,
					p,
					q,
					n,
					w;
				
		bool    	go_right = false;
		
		_SimpleList da ((unsigned long)32);
		
		// try to find the node or where to insert it
		
		for (q=z, p=y; p>=0; q=p, p=go_right?rightChild.lData[p]:leftChild.lData[p])
		{
			long comp = dataList->Compare (b, p);
			if (comp == 0)
			{
				if (cp == false && clear)
					DeleteObject (b);
				return -p-1;
			}
			if (balanceFactor.lData[p] != 0)
			{
				z = q; 
				y = p;
				da.Clear();
			}
			go_right = comp > 0;
			da << go_right;
		}
		
		/*if (da.lLength > 3*log (dataList->lLength+2))
		{
 			WarnError ("AVLList internal error!");
			return -1;
		}*/
		
		// insert new node
		
		n = InsertData (b, xtra,cp);
		
		if (go_right)
			rightChild.lData[q] = n;
		else
			leftChild.lData[q] = n;
			
			
		// update balance factors
			
		p = y;
			
		for (long k=0; p!=n; p=da.lData[k]?rightChild.lData[p]:leftChild.lData[p],k++)
			if (da.lData[k] == 0)
				balanceFactor.lData[p]--;
			else
				balanceFactor.lData[p]++;
				
				
		//if (z < 0)
		//{
			//ConsistencyCheck();
			//return n;
		//}

		if (balanceFactor.lData[y] == -2) 
		{//152
				
			long x = leftChild.lData[y];
			if (balanceFactor.lData[x] == -1) //155
			{
				w 					 = x;
				leftChild.lData [y]  = rightChild.lData[x];
				rightChild.lData[x]  = y;
				balanceFactor.lData[x] = balanceFactor.lData[y] = 0;
			}
			else //156
			{
				w = rightChild.lData[x];
				rightChild.lData[x] = leftChild.lData[w];
				leftChild.lData[w] = x;
				leftChild.lData[y] = rightChild.lData[w];
				rightChild.lData[w] = y;
				if (balanceFactor.lData[w] == -1)
				{
					balanceFactor.lData[x] = 0;
					balanceFactor.lData[y] = 1;
				}
				else
					if (balanceFactor.lData[w] == 0)
					{
						balanceFactor.lData[x] = 0;
						balanceFactor.lData[y] = 0;					
					}
					else
					{
						balanceFactor.lData[x] = -1;
						balanceFactor.lData[y] = 0;
					}
					
				balanceFactor.lData[w] = 0;
			}
		}
		else
			if (balanceFactor.lData[y] == 2)
			{
				long x = rightChild.lData[y];
				if (balanceFactor.lData[x] == 1)
				{
					w 					   = x;
					rightChild.lData [y]   = leftChild.lData[x];
					leftChild.lData[x]     = y;
					balanceFactor.lData[x] = balanceFactor.lData[y] = 0;
				}
				else
				{
					w = leftChild.lData[x];
					leftChild.lData[x] = rightChild.lData[w];
					rightChild.lData[w] = x;
					rightChild.lData[y] = leftChild.lData[w];
					leftChild.lData[w] = y;
					if (balanceFactor.lData[w] == 1)
					{
						balanceFactor.lData[x] = 0;
						balanceFactor.lData[y] = -1;
					}
					else
						if (balanceFactor.lData[w] == 0)
						{
							balanceFactor.lData[x] = 0;
							balanceFactor.lData[y] = 0;					
						}
						else
						{
							balanceFactor.lData[x] = 1;
							balanceFactor.lData[y] = 0;
						}
						
					balanceFactor.lData[w] = 0;
				}				
			}
			else
			{
				//ConsistencyCheck ();
				return n;
			}
			
		if (z >= 0)
		{
			if (y == leftChild.lData[z])
				leftChild.lData[z] = w;
			else
				rightChild.lData[z] = w;
		}
			
		if (y==root)
			root = w;
			
		//ConsistencyCheck ();
						
		return p;
	}
	
	/*dataList->InsertElement (b,-1,false,false);
	leftChild  << -1;
	rightChild << -1;
	balanceFactor << 0;*/
	root =InsertData (b, xtra,cp);
			
	return 0;
		
}

//______________________________________________________________

bool  _AVLList::HasData (long idx)
{
	return leftChild.lData[idx] != 2;
}

//______________________________________________________________

void  _AVLList::Delete (BaseRef b, bool delMe)
{

	if (root == -1)
		return;
		
	_SimpleList pa ((unsigned long)64),
				da ((unsigned long)64);
				
	long		p = root,
				cmp = dataList->Compare (b,p),
				k = 0;
				
	pa.lData[k] = -1;
	da.lData[k++] = 1;
				
	for (;cmp !=0; cmp = dataList->Compare (b,p))
	{
		bool go_right = cmp > 0;
		
		pa.lData[k] = p;
		da.lData[k++] = go_right;
		
		if (go_right)
			p = rightChild.lData[p];
		else
			p = leftChild.lData[p];
			
		if (p<0)
			return;
	}
	
	if (k==1)
		pa.lData[k]   = -1;
	
	emptySlots << p;
	if (delMe)
		DeleteObject (Retrieve(p));
	//((BaseRef*)dataList->lData)[p] = nil;
	dataList->lData[p] = 0;
	DeleteXtra (p);
	
	long r = rightChild.lData[p];

	if (r < 0)
	{
		if (k>1)
		{
			if (da.lData[k-1] == 1)
				rightChild.lData[pa.lData[k-1]] = leftChild.lData[p];
			else
				leftChild.lData[pa.lData[k-1]] = leftChild.lData[p];
		}

		if (p==root)
			//root = pa.lData[k-1];
			root = leftChild.lData[root];
	}
	else
	{
		if (leftChild.lData[r] < 0)
		{
			leftChild.lData[r]     = leftChild.lData[p];
			balanceFactor.lData[r] = balanceFactor.lData[p];
			if (k>1)
			{
				if (da.lData[k-1] == 1)
					rightChild.lData[pa.lData[k-1]] = r;
				else
					leftChild.lData[pa.lData[k-1]] = r;
			}
			else
				root = r;
				
			da.lData[k]   = 1;
			pa.lData[k++] = r;
			//if (p==root)
				//root = r;
		}
		else
		{
			long s;
			int  j = k++;
			for (;;)
			{
				da.lData[k]   = 0;
				pa.lData[k++] = r;
				s = leftChild.lData[r];
				if (leftChild.lData[s] < 0)
					break;
				r = s;
			}
			
			
			leftChild.lData[s] = leftChild.lData[p];
			leftChild.lData[r] = rightChild.lData[s];
			rightChild.lData[s] = rightChild.lData[p];
			balanceFactor.lData[s] = balanceFactor.lData[p];
			
			if (j>1)
			{
				if (da.lData[j-1] == 1)
					rightChild.lData[pa.lData[j-1]] = s;
				else
					leftChild.lData[pa.lData[j-1]] = s;
			}
				
			da.lData[j] = 1;
			pa.lData[j] = s;
			if (p==root)
				root = s;
		}
	}
	
	//if (k>63)
	//{
		//WarnError ("Internal List error");
	//}
	
	while (--k > 0)
	{
		long y = pa.lData[k];
		if (da.lData[k] == 0)
		{
			balanceFactor.lData[y] ++;
			if (balanceFactor.lData[y] == 1)
				break;
			else
				if (balanceFactor.lData[y] == 2)
				{
					long x = rightChild.lData[y];
					if (balanceFactor.lData[x] == -1)
					{
						long w = leftChild.lData[x];
						leftChild.lData[x] = rightChild.lData[w];
						rightChild.lData[w] = x;
						rightChild.lData[y] = leftChild.lData[w];
						leftChild.lData[w] = y;
						if (balanceFactor.lData[w] == 1)
						{
							balanceFactor.lData[x] = 0;
							balanceFactor.lData[y] = -1;
						}
						else
							if (balanceFactor.lData[w] == 0)
							{
								balanceFactor.lData[x] = 0;
								balanceFactor.lData[y] = 0;					
							}
							else
							{
								balanceFactor.lData[x] = 1;
								balanceFactor.lData[y] = 0;
							}
							
						balanceFactor.lData[w] = 0;	
						if (k>1)
						{				
							if (da.lData[k-1] == 1)
								rightChild.lData[pa.lData[k-1]] = w;
							else
								leftChild.lData[pa.lData[k-1]] = w;
						}
						else
							root = w;
							
						//if (y==root)
						//	root = w;
					}
					else
					{
						rightChild.lData[y] = leftChild.lData[x];
						leftChild.lData[x] = y;
						
						if (k>1)
						{
							if (da.lData[k-1] == 1)
								rightChild.lData[pa.lData[k-1]] = x;
							else
								leftChild.lData[pa.lData[k-1]] = x;
						}
						else
							root = x;
							
						if (balanceFactor.lData[x] == 0)
						{
							balanceFactor.lData[x] = -1;
							balanceFactor.lData[y] = 1;
							break;
						}
						else
						{
							balanceFactor.lData[x] = 0;
							balanceFactor.lData[y] = 0;
						}
					}			
				}
		}
		else
		{
			balanceFactor.lData[y] --;
			if (balanceFactor.lData[y] == -1)
				break;
			else	
				if ( balanceFactor.lData[y] == -2)
				{
					long x = leftChild.lData[y];
					if (balanceFactor.lData[x] == 1)
					{
						long w = rightChild.lData[x];
						rightChild.lData[x] = leftChild.lData[w];
						leftChild.lData[w] = x;
						leftChild.lData[y] = rightChild.lData[w];
						rightChild.lData[w] = y;
						if (balanceFactor.lData[w] == -1)
						{
							balanceFactor.lData[x] = 0;
							balanceFactor.lData[y] = 1;
						}
						else
							if (balanceFactor.lData[w] == 0)
							{
								balanceFactor.lData[x] = 0;
								balanceFactor.lData[y] = 0;					
							}
							else
							{
								balanceFactor.lData[x] = -1;
								balanceFactor.lData[y] = 0;
							}
							
						balanceFactor.lData[w] = 0;
						if (k>1)
						{
							if (da.lData[k-1] == 1)
								rightChild.lData[pa.lData[k-1]] = w;
							else
								leftChild.lData[pa.lData[k-1]] = w;
						}
						else
							root = w;
							
						//if (y==root)
							//root = w;
					}
					else
					{
						leftChild.lData[y] = rightChild.lData[x];
						rightChild.lData[x] = y;
						if (k>1)
						{
							if (da.lData[k-1] == 1)
								rightChild.lData[pa.lData[k-1]] = x;
							else
								leftChild.lData[pa.lData[k-1]] = x;
						}
						else
							root = x;
							
						if (balanceFactor.lData[x] == 0)
						{
							balanceFactor.lData[x] = 1;
							balanceFactor.lData[y] = -1;
							break;
						}
						else
						{
							balanceFactor.lData[x] = 0;
							balanceFactor.lData[y] = 0;
						}
					}
				}
		}
	}
	//ConsistencyCheck ();
				
}