/*

HyPhy - Hypothesis Testing Using Phylogenies.

This file defines data and interface classes
used in linking HyPhy with MEGA

Written by SL Kosakovsky Pond; June 2007
Dedicated to Comet (http://www.hyphy.org/comet.jpg)
?/?/1999-06/05/2007

Copyright (C) 1997-2007
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

#ifdef SWIG
%module HyPhy
% {
#include "THyPhy.h"
	%
}
#endif

#ifndef _THYPHY_
#define _THYPHY_

// the progress update/cancel handling function prototype

typedef	 bool _ProgressCancelHandler (char*,int,double);
/* if  the function returns false, then HyPhy will attempt to terminate the task and return
	 HyPhy will also populate three agruments
  -- char * : pointer to an (internally handled by HyPhy) string
			  that describes what task HyPhy is working on at the moment
  -- int    : integer progress value (as in % done)
  -- double : floating point value showing current logL score; could be > 0 if
  			: not applicable (e.g. when doing ancestral state reconstruction)

*/

#define		THYPHY_TYPE_COUNT  3 // how many types are there?

#define		THYPHY_TYPE_STRING 0
#define		THYPHY_TYPE_NUMBER 1
#define		THYPHY_TYPE_MATRIX 2

/* Basic return types supported by the interface object */

class _THyPhyString;
class _THyPhyNumber;
class _THyPhyMatrix;

//_________________________________________________________
class _THyPhyReturnObject // abstract base class
{
public:
	virtual int myType (void) = 0;
	virtual ~_THyPhyReturnObject (void) {};

	// the next three functions are for SWIG to
	// allow the casting of returned base class
	// to a superclass

	_THyPhyString*		castToString	();
	_THyPhyNumber*		castToNumber	();
	_THyPhyMatrix*		castToMatrix	();
};

//_________________________________________________________
class _THyPhyString: public _THyPhyReturnObject // basic string class
{
public:
	_THyPhyString			(const char* = 0, long = 0);
	virtual int myType 		(void) {
		return THYPHY_TYPE_STRING;
	}
	virtual ~_THyPhyString  (void);

	long	  sLength; 		// length of the string
	char*	  sData;		// string characters

};

//_________________________________________________________
class _THyPhyNumber: public _THyPhyReturnObject // basic number class
{
public:
	_THyPhyNumber			(double = 0.0);
	virtual int myType 		(void) {
		return THYPHY_TYPE_NUMBER;
	}
	virtual ~_THyPhyNumber  (void) {}

	double	 nValue;		// the actual number
};

//_________________________________________________________
class _THyPhyMatrix: public _THyPhyReturnObject // basic matrix class
{
public:
	_THyPhyMatrix			(void);
	_THyPhyMatrix			(const long,const long,const double*);

	virtual int myType 		(void) {
		return THYPHY_TYPE_MATRIX;
	}
	virtual ~_THyPhyMatrix   (void);

	double	MatrixCell		(long, long);
	// retrieve the value at a given (0-based!)
	// matrix index

	long	 mRows,
			 mCols; // rows and columns of the matrix

	double*  mData; // the array storing matrix cells
};



//_________________________________________________________
class _THyPhy
{

public:

	_THyPhy	 			 (_ProgressCancelHandler*, const char*, long = 1);
	// create an instance with a given progress/cancel handler
	// the second argument is a path to the directory that
	// HyPhy will use to find its scripting components
	// takes care of needed istantiations (if need)
	// the last argument is relevant for multi-threaded builds
	// and specifies the number of threads that HyPhy should try to spawn


	_THyPhy	 			 (const char*, long = 1);
	// use a default (return true)
	// callback

	~_THyPhy			(void);
	// destructor

	_THyPhyString
	*ExecuteBF			(const char*, bool = true);
	/* take in a string of (null terminated )commands in BF,
	   execute them and return the result (if the BF finished
	   with a return statement)
	   The boolean flag, if true, will clean the previous HyPhy state
	   (suitable for starting a new analysis from scratch)
	   HyPhy will take care of managing the memory needed for the
	   result string
	*/

	void	InitTHyPhy			(_ProgressCancelHandler*, const char*, long);

	void	ClearAll			(void);
	// reset all internal variables, but keep the handler
	// and the file path

	void*	AskFor				(const char*);
	/*
	   After the analysis has finished, this function can be called to
	   retrieve result objects by a string key. Available keys and
	   return mechanisms must be defined via an _THyPhyAskFor function
	   inside the batch commands passed to ExecuteBF.

	   NULL is returned if no result matching the key can be returned,
	   otherwise, the result is a pointer to an object in an internal HyPhy
	   format

	   Memory used by the return object should be de-allocated
	   using DumpResult (see below)

	 */

	void	DumpResult			(void*);

	bool	CanCast				(const void*, const int);
	/* Can the result from AskFor be typecast to the requested type
	   (second argument) - see above for types */

	_THyPhyReturnObject*
	CastResult			(const void*, const int);
	/* Takes a result returned by AskFor, and attempts to represent it
	   in one of the basic return types passed as the 2nd argument
	   (see above).

	   If typecasting fails, NULL is returned */


	void	SetCallbackHandler  (_ProgressCancelHandler*);
	// change the callback handler

	_ProgressCancelHandler *
	GetCallbackHandler (void);
	// return the pointer to the current callback handler


	// for the following 3 functions, the end user is reponsible
	// for disposing of the return structure

	_THyPhyString *  GetWarnings (void) {
		return ConvertHyPhyString (warnings);
	}
	// retrieve warnings (if any) from the last run

	_THyPhyString *  GetErrors	 (void) {
		return ConvertHyPhyString (errors);
	}
	// retrieve errors (if any) from the last run

	_THyPhyString *  GetStdout	 (void) {
		return ConvertHyPhyString (textout);
	}
	// retrieve standard out (of any) from the last run

	// these functions are for internal use only...

	void		PushWarning   (void*);
	void		PushError	  (void*);
	void		PushOutString (void*);

private:

	_THyPhyString			  * ConvertHyPhyString (void*);

	_ProgressCancelHandler	  * theHandler;
	_THyPhyString			  * currentResultHolder,
							  * baseDirectoryInstance;
	// stores the most recent result from ExecuteBF
	long					    askFID;
	// internal ID of _THyPhyAskFor function
	// set by the last call to ExecuteBF
	// -1 is undefined
	void						*errors,
								*warnings,
								*textout;
};

//_________________________________________________________

long	_THyPhyGetLongStatus		(void);
char*	_THyPhyGetStringStatus		(void);
double  _THyPhyGetDoubleStatus		(void);

//void	SetGlobalInterfaceInstance  (const _THyPhy*);

extern	_THyPhy* 					globalInterfaceInstance;
// this object MUST be set to a new instance of _THyPhy
// in order for the rest of the code to intreface with it properly


#endif
