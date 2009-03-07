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

#include <stdlib.h>
#include <stdio.h>
#include "errorfns.h"
#include "hy_strings.h"

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif
	
_String DecodeError (long errCode);

//_____________________________________________________________

int	     gError 				= _HYNOERROR;
bool  	 isFixable 				= true,
		 skipWarningMessages 	= false;
 
//_____________________________________________________________

bool gStatus(void) {
	return gError==_HYNOERROR;
}

//_____________________________________________________________

_String DecodeError (long errCode)
{
	switch (errCode)
	{
		case -101: return "Incompatible Operands";
		break;
		case -102: return "Operation Undefined for Type";
		break;
		case -103: return "Incompatible Matrix Dimensions";
		break;
		case -104: return "Bad Matrix Definition";
		break;
		case -105: return "Matrix Index Out of Range";
		break;
		case -106: return "Bad Matrix Index";
		break;
		case -108: return "Memory Full";
		break;
		case -109: return "Syntax Error";
		break;
		case -110: return "Runtime Expression Error";
		break;
		case -111: return "Non-polynomial expression encountered in polynomial calculation";
		break;
		case -171: return "Dataset index reference out of range";
		break;
		case -200: return "Export Matrix Called With a Non-polynomial Matrix Argument";
		break;
		default:
		return "Unclassified Error";
	}
}

//_____________________________________________________________
bool isError (long errCode) {
//	if (errCode == 0)
//		acknError (false);
	if (errCode!=_HYNOERROR)  {
		if (errCode<0)
		{
			gError = errCode;
			isFixable = TRUE;
			return FALSE;
        }
	}
	gError = _HYNOERROR;
	return TRUE;
}

//_____________________________________________________________

void	warnError (long errCode)
{
	if (errCode == -108)
		warnError (DecodeError (errCode)&_String(" Exiting..."));
	else
		WarnError (DecodeError (errCode)&_String(" Exiting..."));
}

//_____________________________________________________________

void	warnError (const char* theError)
{
	FlagError (theError);
}

//_____________________________________________________________

void	acknError (const char* theError)
{
	WarnError (theError);
}

//_____________________________________________________________

void*	checkPointer (void* p)
{
	if (p)
		return p;
	
	warnError(-108);
	return nil;
}
//_____________________________________________________________
//EOF			
		
	


	
