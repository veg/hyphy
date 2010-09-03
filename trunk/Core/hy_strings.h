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

#ifndef _HSTRINGS_
#define _HSTRINGS_
//#pragma once
#include "baseobj.h"
#include "hy_lists.h"


class _String:public BaseObj {
	
			// contructor/destructor methods
			public:
			
			_String (void); 
				//does nothing
			_String (unsigned long sL, bool flag = false);
				//length constructor
			_String (long);
				//from a long number
			_String (_Parameter);
				//from a floating number
			_String (const _String&);
				// stack copy contructor
			_String (_String*);
				// stack copy contructor
			_String (const _String&, long, long);
				// cut a range of the original string
			_String (const char*);
				// data constructor
			_String (const char);
				// data constructor
			_String (FILE*);
				// read from open file
				
virtual 	~_String(void);
 				//destructor
 				
virtual		BaseRef	makeDynamic (void);
				// create a dynamic copy of this object 				
 				
virtual		void	Initialize (void);

virtual		void	Duplicate (BaseRef);

virtual		void	DuplicateErasing (BaseRef);

const		char	getChar			 (long);
			void	setChar			 (long,char);
			
			void	CopyDynamicString (_String*, bool = true);

virtual 	char& operator [] (long);
  				// element location functions
	 		char operator  () (unsigned long);
 			
 			void operator = (_String);
 				// assignment operator
 			
 			unsigned long Length(void);
 				// string length
 			
 			_String operator & (_String);
  			
 			virtual void operator << (const _String*);
 				// append into operator
	
			void	AppendNewInstance (_String*);
			/* SLKP 20100903: added a utility function to append the contents of a dynamic string 
							  to a string buffer and then free the argument */

 			virtual void operator << (const char);
 				// append into operator
				
			virtual void EscapeAndAppend (const char, char);
				// escape a character and append to this string

			virtual void EscapeAndAppend (const _String &, char mode = 0);
				// escape all characters in a string and append to this string
				// mode = 0 : normal "text" escaping
				// mode = 1: PostScript escaping
				// mode = 2: SQLite escaping

 			virtual void operator << (const char*);
 				// append into operator

 			virtual void Finalize (void);
 				// append into operator finalizer
 				
 			bool ContainsSubstring (_String&);
 				
virtual		BaseRef	toStr (void);

virtual		operator const char* (void);
 				// return good ole char*
 				
 		    char*    getStr(void);
 				
 			_String Chop(long, long);
 				// chop the segment from, to (-1 for any means from beginning/to end)
 				
 			_String Cut (long, long);
 				// cut the string between from and to

			void	Flip(void);
				// in-place flip string
	
			void 	Trim(long, long, bool = false);
 				// trim the string between from and to
 				
 			void	Insert (char, long);
 				// insert a char at a given position (-1 - append)
 				
 			void	Delete (long, long);
 				// delete a range of chars from the string

 			_String Replace(_String, _String, bool);
 				// replace string 1 with string 2, all occurences true/false
 				
 			char FirstNonSpace(long start = 0, long end = -1, char direction = 1);
 				// scans the string for the first non space char

 			long FirstNonSpaceIndex(long start = 0, long end = -1, char direction = 1);
 				// scans the string for the first non space char and returns its index

 			long FirstSpaceIndex(long start = 0, long end = -1, char direction = 1);
 				// scans the string for the first  space char and returns its index

 			long FindEndOfIdent(long start = 0, long end = -1, char wild = '*');
 				// scans the string for the end of an ident starting at 'start'

 			long Find(_String s, long from = 0, long to = -1);
 				// find first occurence of the string between from and to
 				
 			long Find(char s, long from = 0, long to = -1);
 				// find first occurence of the string between from and to

 			long FindAnyCase (_String, long from = 0, long to = -1);
 				// find first occurence of the string between from and to
 				// case insensitive

  			long FindBackwards(_String, long, long);
 				// find first occurence of the string between from and to searching backwards
 				
			long FindBinary(char);
 				// finds a character in an alphabetized string
 				
 			long Adler32 (void);
	
			void FormatTimeString (long);
				// convert the seconds argument into a hh:mm:ss format

 			bool operator == (_String);
 				// lexicographic comparison
 				
 			bool Equal   (_String*);
 				// faster lexicographic comparison

 			char Compare (_String*);
 				// complete lexicographic comparison

 			bool EqualWithWildChar (_String* s, char wildChar = '*');
 				// faster lexicographic comparison

			bool operator > (_String);
 				// lexicographic comparison
 			
 			bool operator < (_String);
 				// lexicographic comparison

			bool Greater (_String*);
 				// lexicographic comparison
 			
 			bool Less (_String*);
 				// lexicographic comparison

			bool operator >= (_String);
 				// lexicographic comparison
 			
 			bool operator <= (_String);
 				// lexicographic comparison

 			bool operator != (_String);
 				// lexicographic comparison
 			
 			bool contains (_String);
 				// superset
 				
			bool contains (char);
 				// superset

 			bool beginswith (_String, bool = true);
 				// begins with string / case sensitive

 			bool startswith (_String&);
 				// begins with string
 				
 			bool endswith (_String, bool = true);
 				// ends with string / case sensitive
 				 				
 			void 	UpCase (void);
 				// upcase the string
 			
 			void 	LoCase (void);
 				// upcase the string
	
			void	AppendAnAssignmentToBuffer (_String*, _String*, bool = true, bool = false, bool = false);
			/* SLKP 20090817
				a utility function to append a statement of the form 
				id = value; to the current string assumed to be in the buffer form
				bool flags have the following meanings:
					3rd: free the 2nd string argument when done
					4th: put quotes around the value
					5th: use := instead of =
 			*/

			void	AppendVariableValueAVL (_String*, _SimpleList&);
			/* SLKP 20090817
				a utility function to append a statement of the form 
				id["varname"] = varvalue; for each variable in the SimpleList arguments
					for String valued variables, their values are properly quoted
			*/

			_List*  Tokenize (_String);
	
  			void    ProcessFileName (bool isWrite = false, bool acceptStringVars = false, Ptr = nil);
 			
  			void    ProcessParameter (void);
  			
  			_String PathComposition (_String);

  			_String PathSubtraction (_String&, char);
			
 			void    StripQuotes(void);
 			
 			bool    IsValidIdentifier(bool = true);
 			
  			bool    IsValidRefIdentifier(void);

			void	ConvertToAnIdent (bool = true);
 			
 			void	KillSpaces		 (_String&);
 			
			void	CompressSpaces	 (void);

 			_String ShortenVarID	 (_String&);

 			static	unsigned long	  storageIncrement;
 			
 			void	RegExpMatch		 (Ptr, _SimpleList&);
 			void	RegExpMatchAll	 (Ptr, _SimpleList&);

  			void	RegExpMatchOnce	 (_String*, _SimpleList&, bool, bool);
  			
  			_String*Sort			 (_SimpleList* = nil);
  			
  			long	LempelZivProductionHistory 
  									 (_SimpleList* = nil);
	
			_Parameter
					ProcessTreeBranchLength ();

	
			/* SLKP 20100831: a utility function to handle 
							  the conversion of branch length strings to parameters
			 
			 */
			 
_Parameter		toNum (void);
 			
			void	SetLength (unsigned long nl) {sLength=nl;}
	
			long	ExtractEnclosedExpression (long&, char, char, bool, bool);
			/* SLKP 20090803
					starting at index [argument 1],
					find a span that encloses an expression (nested) delimited by char[argument 2]
					and char[argument 3] (e.g. {}, ()) respecting quotes (argument 4), and allowing
					escaped characters (argument 5)
			 
					the starting position of the segment will be stored in argument 1 and the
					ending position will be returned. 
			 
					-1 is returned if the starting character could not be found or the expression 
					did not terminate before the end of the string
			*/

			long	FindTerminator			(long, _String&);
			/* SLKP 20090805
			 starting at index [argument 1],
			 find a span that terminates in one of the characters in [argument 2], while
			 respecting (), [], {}, "" and escapes
			 
			 -1 is returned if the starting character could not be found or the expression 
			 did not terminate before the end of the string
			 */
	
 			
 			// data fields
 			unsigned	  long sLength;
 			
 			Ptr			  sData;
 			
};


// _______________________________________________________________________



extern _String empty,
			   emptyAssociativeList,
			   hyphyCiteString;
			   
#ifdef  __MAC__
	extern _String volumeName;
#endif

void	SetStatusBarValue 		    (long,_Parameter,_Parameter);
void	SetStatusLine 			    (_String);

Ptr		PrepRegExp					(_String*, int&, bool);
void	FlushRegExp					(Ptr);
_String GetRegExpError				(int);
void	ReportWarning				(_String);
void	FlagError					(_String);
void	WarnError					(_String);
_String GetVersionString			(void);
_String GetTimeStamp				(bool = false);
void	StringToConsole				(_String&);
void	BufferToConsole				(const char*);
void	NLToConsole					(void);
_String*StringFromConsole			(bool=true);
char	GetPlatformDirectoryChar	(void);

extern  _String 		  __KERNEL__VERSION__;


#endif
 			