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

*/

#include <math.h>
#include <float.h>
#include <limits.h>

#include "likefunc.h"
#include "parser.h"
#include "matrix.h"
#include "stdlib.h"
#include "string.h"
#ifdef   __MACHACKMP__
	#undef	 __MACHACKMP__
	#define  _BSD_CLOCK_T_DEFINED_
	#include <stdio.h>
	#define  __MACHACKMP__
	#include "time.h"
#else
	#include <stdio.h>
	#include "time.h"
#endif
#include "ctype.h"
#include "polynoml.h"
#include "batchlan.h"

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif

_Formula		*chi2 = nil, 
				*derchi2 = nil;

				
extern	
_SimpleList		simpleOperationCodes, 
				simpleOperationFunctions;

_List			globalNamesSupportList,
				varNamesSupportList,
		 		variablePtrs;   // stores all the variables declared so far
		 		
_AVLList		*lookAside = nil;
		 
_AVLListX   	variableNames (&varNamesSupportList),
				_hyApplicationGlobals (&globalNamesSupportList);

_Parameter 	  	sqrtPi 						= 1.77245385090551603, 
				maxRombergSteps 			= 8., 
			 	integrationPrecisionFactor	= 1.e-5,
			 	printDigits;
			 	
_String			intPrecFact ("INTEGRATION_PRECISION_FACTOR"),
				intMaxIter  ("INTEGRATION_MAX_ITERATES");

long			lastMatrixDeclared = -1, 
				dummyVariable1, 
				dummyVariable2,
				expressionsParsed = 0;
	
// indices of all independent variables
void 			DeleteVariable (long dv);

_List			FunctionNameList,
				BuiltInFunctions;
				
_SimpleList		FunctionArgumentCount,
				BuiltInFunctionParameterCount,
				freeSlots,
			  * deferSetFormula = nil,
				deferIsConstant;
				
bool			useGlobalUpdateFlag = false;

#ifndef  __HYALTIVEC__
	_Parameter  machineEps = 1e-12,
				tolerance  = DBL_EPSILON;
#else
	_Parameter  machineEps = 1e-7,
				tolerance  = FLT_EPSILON;
#endif

_Parameter	 	gammaCoeff [7] = {
					 2.50662827463100050,
					 190.9551718944012,
					 -216.8366818451899,
					 60.19441758801798,
					 -3.087513097785903,
					 0.003029460875352382,
					 -0.00001345152485367085}, 
					 
				lngammaCoeff [6] = {
					76.18009172947146,
					-86.50532032941677,
					24.01409824083091,
					-1.231739572450155,
					0.1208650973866179e-2,
					-0.5395239384953e-5},
				
				pi_const = 3.141592653589793,
				long_max = (_Parameter)LONG_MAX;

//__________________________________________________________________________________
BaseRef			parameterToString 		(_Parameter);
void 			DeleteTreeVariable 		(long, _SimpleList &, bool);


//__________________________________________________________________________________
_Variable * LocateVar (long index)
{
	return (_Variable *)(((BaseRef*)variablePtrs.lData)[index]);
}

//__________________________________________________________________________________
BaseRef		parameterToString (_Parameter value)
{
	char dump [255];
	long digs = printDigits;
	if (digs<=0 || digs>15)
	{
		if (round(value) == value && fabs (value) < long_max)
			sprintf (dump,"%ld",lrint (value));
		else
			sprintf (dump,PRINTF_FORMAT_STRING,value);
	}
	else
	{
		_String format("%-");
		#ifdef __USE_LONG_DOUBLE__
			format = format&_String(digs)&"Lg";
		#else		
			format = format&_String(digs)&'g';
		#endif
		sprintf (dump,(const char*)format.sData,value);
	}		
	return new _String (dump);
}

//__________________________________________________________________________________
_PMathObj	FetchObjectFromVariableByType (_String* id, int objectClass)
{
	if (id)
	{
		_Variable * v = FetchVar (LocateVarByName (*id));
		if (v && (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass))
			return v->Compute();
	}
	return nil;
}

//__________________________________________________________________________________
_PMathObj	FetchObjectFromVariableByTypeIndex (long idx, int objectClass)
{
	_Variable * v = FetchVar (idx);
	if (v && (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass))
		return v->GetValue();
	return nil;
}

//__________________________________________________________________________________
long LocateVarByName (_String& name)
{
	return variableNames.Find (&name);
}


//__________________________________________________________________________________
_Variable* FetchVar (long index)
{
	return index>=0?(_Variable *)variablePtrs(variableNames.GetXtra(index)):nil;
}

//__________________________________________________________________________________
void	   UpdateChangingFlas (long vN)
{
	// check to see if formulae contain a reference to this var
	// if so: "decompile" them
	
	long topLimit = compiledFormulaeParameters.lLength;
		
	_SimpleList * toDelete = nil;
	
	for (long k = 0; k<topLimit; k++)
	{
		long g = ((_SimpleList*)compiledFormulaeParameters.lData[k])->BinaryFind (vN,0);
		
		if (g>=0)
		{
			_ElementaryCommand* thisCommand = (_ElementaryCommand*)listOfCompiledFormulae.lData[k];
			_Formula  *f  = (_Formula*)(thisCommand->simpleParameters.lData[1]),
				  	  *f2 = (_Formula*)(thisCommand->simpleParameters.lData[2]);
				  	  
			delete f;
			delete f2;
			
			thisCommand->simpleParameters.Clear();
			
			// MOD 10/21/2005
			
			if (!toDelete)
				checkPointer(toDelete = new _SimpleList);
			
			//listOfCompiledFormulae.    Delete(k);
			//compiledFormulaeParameters.Delete(k);
			
			//k--;
			//topLimit--;
			
			*toDelete << k;
		}
	}
	
	if (toDelete)
	{
		listOfCompiledFormulae.DeleteList 		(*toDelete);
		compiledFormulaeParameters.DeleteList	(*toDelete);
		DeleteObject (toDelete);
	}
}

//__________________________________________________________________________________
void	   UpdateChangingFlas (_SimpleList & involvedVariables)
{
	
	long 		  topLimit 		   = compiledFormulaeParameters.lLength;
	_SimpleList * toDelete 		   = nil;
	
	for (long k = 0; k<topLimit; k++)
	{
		long g = ((_SimpleList*)compiledFormulaeParameters.lData[k])->CountCommonElements (involvedVariables,true);
		
		if (g>0)
		{
			_ElementaryCommand* thisCommand = (_ElementaryCommand*)listOfCompiledFormulae.lData[k];
			
			_Formula  *f  = (_Formula*)(thisCommand->simpleParameters.lData[1]),
				  	  *f2 = (_Formula*)(thisCommand->simpleParameters.lData[2]);
				  	  
			delete f;
			delete f2;
			
			thisCommand->simpleParameters.Clear();
			
			if (!toDelete)
				checkPointer(toDelete = new _SimpleList);
			
			*toDelete << k;
		}
	}
	
	if (toDelete)
	{
		listOfCompiledFormulae.DeleteList 		(*toDelete);
		compiledFormulaeParameters.DeleteList	(*toDelete);
		DeleteObject (toDelete);
	}
}

//__________________________________________________________________________________
void DeleteVariable (long dv, bool deleteself)
{
	if (dv>=0)
	{	
	
		_String *name  = (_String*)variableNames.Retrieve (dv);
		_String myName = *name&'.';
		long	vidx   = variableNames.GetXtra (dv);
		
		UpdateChangingFlas (vidx);
		
		_SimpleList recCache;
		 variableNames.Find (name,recCache);
		_String    	nextVarID;// = *(_String*)variableNames.Retrieve(variableNames.Next (dv,recCache));
		long	    nvid;
		if ((nvid = variableNames.Next (dv,recCache))>=0)
			nextVarID = *(_String*)variableNames.Retrieve(nvid);
		
		if (deleteself)
		{	
			_SimpleList tcache;
			long		iv,
						k = variableNames.Traverser (tcache, iv, variableNames.GetRoot());
						
			for (; k>=0; k = variableNames.Traverser (tcache, iv))
			{
				_Variable * thisVar = FetchVar(k);
				
				if (thisVar->CheckFForDependence (vidx,false))
				{	
					_PMathObj curValue = thisVar->Compute();
					curValue->nInstances++; // this could be a leak 01/05/2004.
					thisVar->SetValue (curValue);
					DeleteObject (curValue);
				}
			}		

			_Variable* delvar = (FetchVar(dv));
			DeleteObject (delvar);
			
			variableNames.Delete (variableNames.Retrieve(dv),true);
			(*((_SimpleList*)&variablePtrs))[vidx]=0;
			freeSlots<<vidx;
		}
		else
		{
			_Variable* delvar = (FetchVar(dv));
			if (delvar->IsContainer())
			{
				_VariableContainer* dc = (_VariableContainer*)delvar;
				dc->Clear();
			}
		}
		
		_List		toDelete;
					
		recCache.Clear();
		long nextVar = variableNames.Find (&nextVarID,recCache);

		for (; nextVar>=0; nextVar = variableNames.Next (nextVar, recCache))
		{
			_String dependent = *(_String*)variableNames.Retrieve (nextVar);
			if (dependent.startswith(myName))
				toDelete && & dependent;
			else
				break;
		}
		
		for (long k=0; k< toDelete.lLength; k++)
			DeleteVariable (*(_String*)toDelete(k));
	}	
}

//__________________________________________________________________________________
void DeleteTreeVariable (long dv, _SimpleList & parms, bool doDeps)
{
	if (dv>=0)
	{	
		_String *name  = (_String*)variableNames.Retrieve (dv);
		_String myName = *name&".";
		long	vidx   = variableNames.GetXtra (dv);

		UpdateChangingFlas (vidx);
		
		_SimpleList recCache;
		variableNames.Find (name,recCache);
		_String    	nextVarID;
		long	    nvid;
		if ((nvid = variableNames.Next (dv,recCache))>=0)
			nextVarID = *(_String*)variableNames.Retrieve(nvid);


		{
			_SimpleList tcache;
			long		iv,
						k = variableNames.Traverser (tcache, iv, variableNames.GetRoot());
						
			for (; k>=0; k = variableNames.Traverser (tcache, iv))
			{
				_Variable * thisVar = FetchVar(k);
				
				if (thisVar->CheckFForDependence (vidx,false))
				{	
					_PMathObj curValue = thisVar->Compute();
					curValue->nInstances++;
					thisVar->SetValue (curValue);
					DeleteObject (curValue);
				}
			}	
		}		

		_Variable* delvar = (FetchVar(dv));
		if (delvar->ObjectClass() != TREE)
		{
			variableNames.Delete (variableNames.Retrieve(dv),true);
			(*((_SimpleList*)&variablePtrs))[vidx]=0;
			freeSlots<<vidx;
			DeleteObject (delvar);
		}
		else
		{	
			((_VariableContainer*)delvar)->Clear();
		}
		if (doDeps)
		{
			_List toDelete;
			recCache.Clear();
			long nextVar = variableNames.Find (&nextVarID,recCache);
			for (; nextVar>=0; nextVar = variableNames.Next (nextVar, recCache))
			{
				_String dependent = *(_String*)variableNames.Retrieve (nextVar);
				if (dependent.startswith(myName))
				{
					if (dependent.Find ('.', myName.sLength+1, -1)>=0)
					{
						_Variable * checkDep = FetchVar (nextVar);
						if (!checkDep->IsIndependent())
						{
							_PMathObj curValue = checkDep->Compute();
							curValue->nInstances++;
							checkDep->SetValue (curValue);
							DeleteObject (curValue);
						}
						parms << variableNames.GetXtra (nextVar);
					}
					else
						toDelete && & dependent;
				}
				else
					break;
			}
			
			for (long k=0; k<toDelete.lLength; k++)
			{
				//StringToConsole (*(_String*)toDelete(k));
				//BufferToConsole ("\n");
				DeleteTreeVariable (*(_String*)toDelete(k),parms,false);
			}
		}
	}	
}
//__________________________________________________________________________________
void DeleteVariable (_String&name, bool deleteself)
{
	DeleteVariable(LocateVarByName (name), deleteself);
}

//__________________________________________________________________________________
void DeleteTreeVariable (_String&name, _SimpleList& parms, bool doDeps)
{
	DeleteTreeVariable(LocateVarByName (name), parms,doDeps);
}

//__________________________________________________________________________________
_Variable* CheckReceptacle (_String* name, _String fID, bool checkValid)
{
	if (checkValid && (!name->IsValidIdentifier()))
	{
		_String errMsg = *name & " is not a variable identifier in call to " & fID;
		WarnError (errMsg);
		return nil;			
	}

	long 	f = LocateVarByName (*name);
	if (f<0)
	{
		_Variable dummy (*name);
		f = LocateVarByName (*name);
	}
	
	return FetchVar(f);
}

//__________________________________________________________________________________
bool CheckReceptacleAndStore (_String* name, _String fID, bool checkValid, _PMathObj v, bool dup)
{
	_Variable * theV = CheckReceptacle(name, fID, checkValid);
	if (theV)
	{
		theV->SetValue (v, dup);
		return true;
	}
	if (!dup)
		DeleteObject (v);
	return false;
}


//__________________________________________________________________________________
void  InsertVar (_Variable* theV)
{
	long pos = variableNames.Insert (theV->theName);
	
	/*if (theV->GetName()->Equal (&_String("PS_2")))
	{
		printf ("Making...\n");
	}*/
	
	if (pos < 0 && isDefiningATree > 1) 
	// automatically fix duplicate autogenerated tree node name
	{
		long trySuffix  = 1;
		_String * tryName = new _String;
		do
		{
			*tryName = *theV->theName & "_" & trySuffix;
			pos      = variableNames.Insert (tryName);
			trySuffix ++;
		}
		while (pos < 0);
		DeleteObject(theV->theName);
		theV->theName = tryName;
	}
	
	if (pos < 0)
	{
		if (isDefiningATree == 1)
		{
			_String errMsg (*theV->GetName());
			errMsg = errMsg& " is already being used - please rename one of the two variables.";
			WarnError(errMsg);
		}
			
		theV->theIndex = variableNames.GetXtra(-pos-1);		
		return;
	}
	else
		theV->theName->nInstances++;

	if (freeSlots.lLength)
	{
		theV->theIndex = freeSlots.lData[freeSlots.lLength-1];
		variablePtrs[theV->theIndex]=theV->makeDynamic();
		freeSlots.Delete(freeSlots.lLength-1);
	}
	else
	{
		theV->theIndex = variablePtrs.lLength;
		variablePtrs&&theV;
	}
	variableNames.SetXtra (pos, theV->theIndex);
}

//__________________________________________________________________________________
_String&  AppendContainerName (_String& inString, _VariableContainer* theP)
{
	static _String returnMe;
	
	if (!theP)
		return inString;
	
	returnMe = *theP->GetName() & '.' & inString;
	return returnMe;
}

//__________________________________________________________________________________
void  RenameVariable (_String* oldName, _String* newName)
{
	_String		oldNamePrefix (*oldName&'.'),
				newNamePrefix (*newName&'.');

	_List 			toRename;
	_SimpleList		xtras,
					traverser;
	
	long f = variableNames.Find (oldName, traverser);		
	if (f>=0)
	{
		toRename << oldName;
		xtras    << variableNames.GetXtra (f);
		f = variableNames.Next (f, traverser);
						
		for	 (;f>=0 && ((_String*)variableNames.Retrieve (f))->startswith (oldNamePrefix); f = variableNames.Next (f, traverser))
		{
			toRename << variableNames.Retrieve (f);
			xtras << variableNames.GetXtra (f);
		}
	}
	
	for (long k = 0; k < toRename.lLength; k++)
	{
		_Variable * thisVar = FetchVar (xtras.lData[k]);
		thisVar->GetName()->nInstances --;
		if (k)
			thisVar->theName = new _String(thisVar->GetName()->Replace(oldNamePrefix,newNamePrefix,true));
		else
			thisVar->theName = new _String(*newName);
			
		variableNames.Delete (toRename (k), true);
		variableNames.Insert (thisVar->GetName(),xtras.lData[k]);
		thisVar->GetName()->nInstances++;
	}
}

//__________________________________________________________________________________
void  ReplaceVar (_Variable* theV)
{
	long pos = variableNames.Find (theV->theName);	
	if (pos>=0)
	{
		pos = variableNames.GetXtra(pos);				
		UpdateChangingFlas 	 (pos);
		variablePtrs.Replace (pos,theV,true);
	}
	else
		InsertVar (theV);	
}

//__________________________________________________________________________________

bool _MathObject::IsDefined	(_String& s)  // is this operation defined for the type
{
	return s.sLength;
}

//__________________________________________________________________________________

_String		UnOps ("-,!,Abs,Sin,Cos,Tan,Exp,Log,Arctan,Time,Gamma,Transpose,Sqrt,Erf,Rows,Columns,LUDecompose,Inverse,BranchCount,TipCount,ZCDF,Eigensystem,Simplex,Type,"), 
			HalfOps (":<>=!&|");
		
_SimpleList	opPrecedence, 
			BinOps,
			associativeOps;
			
//__________________________________________________________________________________

void	SetupOperationLists (void)
{
	BinOps<<'|'*256+'|';
	opPrecedence<<1;
	BinOps<<'&'*256+'&';
	opPrecedence<<2;
	BinOps<<'='*256+'=';
	opPrecedence<<3;
	BinOps<<'!'*256+'=';
	opPrecedence<<3;
	BinOps<<'<';
	opPrecedence<<4;
	BinOps<<'>';
	opPrecedence<<4;
	BinOps<<'<'*256+'=';
	opPrecedence<<4;
	BinOps<<'>'*256+'=';
	opPrecedence<<4;
	BinOps<<'+';
	associativeOps << opPrecedence.lLength;
	opPrecedence<<5;
	BinOps<<'-';
	opPrecedence<<5;
	BinOps<<'*';
	associativeOps << opPrecedence.lLength;
	opPrecedence<<6;
	BinOps<<'/';
	opPrecedence<<6;
	BinOps<<'%';
	opPrecedence<<6;
	BinOps<<'$';
	opPrecedence<<6;
	BinOps<<'^';
	opPrecedence<<7;
	
	_String fName ("Beta");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("BranchLength");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("BranchName");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("CChi2");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("CGammaDist");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 3;
	fName = ("Format");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 3;
	fName = ("GammaDist");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 3;
	fName = ("IBeta");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 3;
	fName = ("IGamma");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("InvChi2");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("LUSolve");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("Max");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("Min");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("PSTreeString");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 3;
	fName = ("Random");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("RerootTree");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("TEXTreeString");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = ("TipName");
	FunctionNameList&& &fName;
	FunctionArgumentCount<< 2;
	fName = "_x_";
	_Variable d(fName,true);
	dummyVariable1 =  variableNames.GetXtra(LocateVarByName (fName));
	fName = "_n_";
	_Variable d1 (fName,true);
	dummyVariable2 =  variableNames.GetXtra(LocateVarByName (fName));
	

	if (BuiltInFunctions.lLength==0)
	// construct a list of operations
	// don't forget to update SimplifyConstants, simpleOperationCodes, InternalDifferentiate, InternalSimplify, Formula::HasChanged and all Execute commands
	// also MAccess and MCoord codes are used in Parse to merge multiple matrix access operations
	{
		_String s;
		s = "!";//0
		BuiltInFunctions&& &s;
		BuiltInFunctionParameterCount<<1;		
		s = "!=";//1
		BuiltInFunctions&& &s;
		s = "$";//2
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "%";//3
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "&&";//4
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "*";//5
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<5;
		simpleOperationFunctions<<(long)MultNumbers;
		s = "+";//6
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<6;
		simpleOperationFunctions<<(long)AddNumbers;
		s = "-";//7
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<7;
		simpleOperationFunctions<<(long)SubNumbers;
		s = "/";//8
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<8;
		simpleOperationFunctions<<(long)DivNumbers;
		s = "<";//9
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<9;
		simpleOperationFunctions<<(long)LessThan;
		s = "<=";//10
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<10;
		simpleOperationFunctions<<(long)LessThanE;
		s = "==";//11
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<11;
		simpleOperationFunctions<<(long)EqualNumbers;
		s = ">";//12
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<12;
		simpleOperationFunctions<<(long)GreaterThan;
		s = ">=";//13
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<13;
		simpleOperationFunctions<<(long)GreaterThanE;
		s = "Abs";//14
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Arctan";//15
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Beta";//16
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "BranchCount";//17
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "BranchLength";//18
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "BranchName";//19
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "CChi2";//20
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "CGammaDist";//21
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "Columns";//22
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Cos";//23
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Eigensystem";//24
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Erf";//25
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Exp";//26
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<26;
		simpleOperationFunctions<<(long)ExpNumbers;
		s = "Format";//27
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "Gamma";//28
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "GammaDist";//29
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "IBeta";//30
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "IGamma";//31
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "InvChi2";//32
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Inverse";//33
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "LUDecompose";//34
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "LUSolve";//35
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "Log";//36
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<36;
		simpleOperationFunctions<<(long)LogNumbers;
		s = "MAccess";//37
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<37;
		simpleOperationFunctions<<(long)FastMxAccess;
		s = "MCoord";//38
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "Max";//39
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<39;
		simpleOperationFunctions<<(long)MaxNumbers;
		s = "Min";//40
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<40;
		simpleOperationFunctions<<(long)MinNumbers;
		s = "PSTreeString";//41
		BuiltInFunctionParameterCount<<3;		
		BuiltInFunctions&& &s;
		s = "Random";//42
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		simpleOperationCodes<<42;
		simpleOperationFunctions<<(long)RandomNumber;
		s = "RerootTree";//43
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "Rows";//44
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Simplex";//45
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Sin";//46
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Sqrt";//47
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "TEXTreeString";//48
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Tan";//49
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Time";//50
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "TipCount";//51
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "TipName";//52
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Transpose";//53
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "Type";//54
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "ZCDF";//55
		BuiltInFunctionParameterCount<<1;		
		BuiltInFunctions&& &s;
		s = "^";//56
		simpleOperationCodes<<56;
		simpleOperationFunctions<<(long)Power;
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
		s = "||";//57
		BuiltInFunctionParameterCount<<2;		
		BuiltInFunctions&& &s;
	}			
}
//__________________________________________________________________________________


_PMathObj _MathObject::Type (void)  
{
	_FString * ts = new _FString();
	switch (ObjectClass())
	{

		case NUMBER:
			*(ts->theString)="Number"; break;
		case MATRIX:
			*(ts->theString)="Matrix"; break;
		case CONTAINER:
			*(ts->theString)="Container"; break;
		case TREE_NODE:
			*(ts->theString)="TreeNode"; break;
		case TREE:
			*(ts->theString)="Tree"; break;
		case STRING:
			*(ts->theString)="String"; break;
		case ASSOCIATIVE_LIST:
			*(ts->theString)="AssociativeList"; break;
		case TOPOLOGY:
			*(ts->theString)="Topology"; break;
		case POLYNOMIAL:
			*(ts->theString)="Polynomial"; break;
		default:
			*(ts->theString) = "Unknown";
				
	}
	
	return ts;
}


//__________________________________________________________________________________
	
	
_PMathObj _MathObject::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{
	switch (opCode)
	{
		case 0: // !
			return LNot();
			break;
		case 1: // !=
			return NotEqual(p);
			break;
		case 2: // $
			return longDiv(p);
			break;
		case 3: // %
			return lDiv(p);
			break;
		case 4: // &&
			return LAnd(p);
			break;
		case 5: // *
			return Mult(p);
			break;
		case 6: // +
			return Add(p);
			break;
		case 7: // -
			if (p)
				return Sub(p);
			else
				return Minus();
			break;
		case 8: // /
			return Div(p);
			break;
		case 9: // <
			return Less(p);
			break;
		case 10: // <=
			return LessEq(p);
			break;
		case 11: // ==
			return AreEqual(p);
			break;
		case 12: // >
			return Greater(p);
			break;
		case 13: // >=
			return GreaterEq(p);
			break;
		case 14: // Abs
			return Abs();
			break;
		case 15: // Arctan
			return Arctan();
			break;
		case 16: // Beta
			return Beta(p);
			break;
		case 20: // CChi2
			return CChi2(p);
			break;
		case 21: // CGammaDist
			return CGammaDist(p,p2);
			break;
		case 23: // Cos
			return Cos();
			break;
		case 22: // Columns
		case 44: // Roes
			return new _Constant (0.0);
			break;
		case 25: // Erf
			return Erf();
			break;
		case 26: // Exp
			return Exp();
			break;
		case 27: // Format
			return FormatNumberString(p,p2);
			break;
		case 28: // Gamma
			return Gamma();
			break;
		case 29: // GammaDist
			return GammaDist(p,p2);
			break;
		case 30: // IBeta
			return IBeta(p,p2);
			break;
		case 31: // IGamma
			return IGamma(p);
			break;
		case 32: // InvChi2
			return InvChi2(p);
			break;
		case 36: // Log
			return Log();
			break;
		case 39: // Max
			return Max(p);
			break;
		case 40: // Min
			return Min(p);
			break;
		case 42: // Random
			return Random(p);
			break;
		case 45: // Simplex
			return Simplex();
			break;
		case 46: // Sin
			return Sin();
			break;
		case 47: // Sqrt
			return Sqrt();
			break;
		case 49: // Tan
			return Tan();
			break;
		case 50: // Time
			return Time();
			break;
		case 54: // Type
			return Type();
			break;
		case 55: // ZCDF
			return ZCDF();
			break;
		case 56: // ^
			return Raise(p);
			break;
		case 57: // ||
			return LOr(p);
			break;
		default:
		{
			_String *ss,errMsg;
			ss = (_String*)toStr();
			errMsg = _String("Operation ")&*(_String*)BuiltInFunctions(opCode)&" is not defined for "&*ss;
			DeleteObject(ss);
			WarnError (errMsg);
			return new _Constant (0.0);
		}
	}
	return (_PMathObj)makeDynamic();
}
	
	
			

//__________________________________________________________________________________

_Stack::_Stack (void)
{
}

//__________________________________________________________________________________

void	_Stack::Initialize (void)
{
	theStack.Initialize();
}

//__________________________________________________________________________________

void	_Stack::Duplicate (BaseRef s)
{
	theStack.Duplicate(&((_Stack*)s)->theStack);
}


//__________________________________________________________________________________

_Stack::~_Stack (void)
{
}

//__________________________________________________________________________________

bool	 _Stack::Push (_PMathObj newObj) 	// push object onto the stack
{
	theStack<<(newObj);
	return true;
}

//__________________________________________________________________________________

_PMathObj _Stack::Pop (bool del)		// pop object from the top of the stack
{
	
	_PMathObj r = (_PMathObj)theStack.lData[theStack.lLength-1];
	if (del)
		theStack.lLength--;
	
	return r;
}

//__________________________________________________________________________________

long	 _Stack::StackDepth (void)	// returns the depth of the stack
{
	return theStack.lLength;
}

//__________________________________________________________________________________

void	 _Stack::Reset (void)	// clears the stack
{
	theStack.Clear();
}
	
	
//__________________________________________________________________________________

_Constant::_Constant (_Parameter value)
{
	theValue = value;
}
//__________________________________________________________________________________

void _Constant::Initialize (void)
{
	BaseObj::Initialize();
	theValue = 0;
}
//__________________________________________________________________________________

void _Constant::Duplicate (BaseRef c)
{
	BaseObj::Initialize();
	theValue = ((_Constant*)c)->theValue;
}

//__________________________________________________________________________________

BaseRef _Constant::makeDynamic (void)
{
	_Constant * res = (_Constant*)checkPointer(new _Constant);
	res->Duplicate(this);
	return res;
}

//__________________________________________________________________________________

_Constant::_Constant (_String& s)
{
	theValue = atof (s.sData);
}

//__________________________________________________________________________________
_Constant::_Constant (void)  {
	theValue = 0;
}

//__________________________________________________________________________________
//_Constant::~_Constant (void)  {
//}
	
//__________________________________________________________________________________
_Parameter	  _Constant::Value (void)
{
	return theValue;
}
//__________________________________________________________________________________
BaseRef _Constant::toStr(void)
{
	return parameterToString(Value());
}
	
//__________________________________________________________________________________
_PMathObj _Constant::Add (_PMathObj theObj) 
{		
	if (theObj->ObjectClass() == STRING)
		return new _Constant ((theValue+((_FString*)theObj)->theString->toNum()));
	else
		return new _Constant ((theValue+((_Constant*)theObj)->theValue));		
}

//__________________________________________________________________________________
_PMathObj _Constant::Sub (_PMathObj theObj) 
{
	//if (theObj) return nil;
	return new _Constant ((theValue-((_Constant*)theObj)->theValue));
	//else
	//	return  nil;
	//return 	   (_PMathObj)result.makeDynamic();
}	

//__________________________________________________________________________________
_PMathObj _Constant::Minus (void) 
{
	return 	   new 	_Constant (-Value());
}
//__________________________________________________________________________________
_PMathObj _Constant::Mult (_PMathObj theObj) 
{
//	if (!theObj) return nil;
	return new _Constant ((theValue*((_Constant*)theObj)->theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Div (_PMathObj theObj) 
{
//	if (!theObj) return nil;
	return new _Constant ((theValue/((_Constant*)theObj)->theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::lDiv (_PMathObj theObj) // %
{
	if (theObj) 
	{
		long 	   denom = ((_Constant*)theObj)->theValue;
		return 	   denom?new _Constant  ((long)(Value())%denom):new _Constant  ((long)(Value()));
	}
	else	
		return nil;	
}
//__________________________________________________________________________________
_PMathObj _Constant::longDiv (_PMathObj theObj) // div
{
	if (theObj) 
	{
		long 	   denom = ((_Constant*)theObj)->theValue;
		return 	   denom?new _Constant  ((long)(Value())/denom):new _Constant  (0.0);
	}
	else	
		return nil;	
}
//__________________________________________________________________________________
_PMathObj _Constant::Raise (_PMathObj theObj) 
{
	if (!theObj) 
		return nil;
		
	_Parameter    base  = Value(),
				  expon = theObj->Value();
	
	if (base>0.0)
	{
		return 	  new  _Constant (exp (log(base)*(expon)));;
	}
	else
	{
		if (base<0.0)
			if (CheckEqual (expon, (long)expon))
				return new _Constant (((((long)expon)%2)?-1:1)*exp (log(-base)*(expon)));
			else
			{
				_String errMsg ("An invalid base/exponent pair passed to ^");
				WarnError (errMsg.sData);
			}
				
		return 	   new _Constant (0.0);
	}
}

//__________________________________________________________________________________

long randomCount = 0;

//__________________________________________________________________________________
_PMathObj _Constant::Random (_PMathObj upperB)
{
	if (randomCount == 0) 
	{
		randomCount++;
	}
	_Parameter l = theValue, u=((_Constant*)upperB)->theValue,r = l;
	if (u>l)
	{
		r=genrand_int32();
		r/=RAND_MAX_32;
		r =l+(u-l)*r;
	}
	return new _Constant (r);
	
}

//__________________________________________________________________________________
void	 _Constant::Assign (_PMathObj theObj) 
{
	this->~_Constant ();
	theValue = ((_Constant*)theObj)->theValue;
}

//__________________________________________________________________________________
bool	 _Constant::Equal (_PMathObj theObj) 
{
	return theValue==((_Constant*)theObj)->theValue;
}

//__________________________________________________________________________________
_PMathObj _Constant::Abs (void) 
{
	return 	   new _Constant (fabs(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Sin (void) 
{
	return 	   new 	_Constant (sin(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Cos (void) 
{
	return 	   new _Constant  (cos(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Tan (void) 
{
	return 	   new _Constant  (tan(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Exp (void) 
{
	return 	   new _Constant  (exp(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::FormatNumberString (_PMathObj p, _PMathObj p2) 
{
	long 	   a1 = p->Value(), 
			   a2 = p2->Value();
			   
	char	   format[32],
			   buffer[255];
			   
	#ifdef	   __USE_LONG_DOUBLE__
	if (a1>=0 && a2>=0)
	{
		if (a1>0)
			sprintf	   (format,"%%%ld.%ldLf",(long)a1,(long)a2);
		else
			sprintf	   (format,"%%.%ldLf",(long)a2);
	}
	else
		if (a1>=0)
			sprintf	   (format,"%%%ldLf",(long)a1);
		else
			if (a2>=0)
				sprintf	   (format,"%%.%ldLf",(long)a2);
			else
				sprintf	   (format,"%%Lg");
	#else
	if (a1>=0 && a2>=0)
	{
		if (a1>0)
			sprintf	   (format,"%%%ld.%ldf",(long)a1,(long)a2);
		else
			sprintf	   (format,"%%.%ldf",(long)a2);
	}
	else
		if (a1>=0)
			sprintf	   (format,"%%%ldf",(long)a1);
		else
			if (a2>=0)
				sprintf	   (format,"%%.%ldf",(long)a2);
			else
				sprintf	   (format,"%%g");
	
	#endif
	a1 = sprintf    (buffer,format,Value());
	_String    t (buffer);
	return 	   new _FString (t);
}
//__________________________________________________________________________________
_PMathObj _Constant::Log (void) 
{
	return 	   new _Constant  (log(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Sqrt (void) 
{
	return 	   new _Constant  (sqrt(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Arctan (void) 
{
	return 	   new _Constant  (atan(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Gamma (void) 
{
	_Parameter theV = theValue>=1.0?theValue:2-theValue, result = gammaCoeff[0], temp = theV;
	
	for (long i=1; i<7; i++, temp+=1.0)
	{
		result+=gammaCoeff[i]/temp;
	}
	
	temp = theV+4.5;
	result *= exp(-temp+log(temp)*(theV-.5));
	
	if (theValue>=1.0)
	{
		return 	  new _Constant  (result);
	}
	
	else
	{
		temp = pi_const*(1-theValue);
		
		return 	   new _Constant  (temp/result/sin(temp));
	}
	return nil;
}

//__________________________________________________________________________________
_PMathObj _Constant::LnGamma (void)
{
	// obtained from Numerical Recipes in C, p. 214 by afyp, February 7, 2007
	_Parameter	x, y, tmp, ser;
	
	y = x = theValue;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	
	for (long j = 0; j <= 5; j++) ser += lngammaCoeff[j] / ++y;
	
	return new _Constant (-tmp + log(2.506628274631005*ser/x));
}

//__________________________________________________________________________________
_PMathObj _Constant::Beta (_PMathObj arg) 
{
	if (arg->ObjectClass()!=NUMBER)
	{
		_String errMsg ("A non-numerical argument passed to Beta(x,y)");
		WarnError (errMsg.sData);
		return	  nil;
	}
	_Constant argVal   = ((_Constant*)arg)->theValue;
	_Constant *result  = (_Constant *)Gamma(), 
			  *result1 = (_Constant *)argVal.Gamma();
			  
	argVal.SetValue (theValue+argVal.theValue);
	_Constant *result2 = (_Constant *)argVal.Gamma();
	argVal.SetValue (result->theValue*result1->theValue/result2->theValue);
	DeleteObject (result);
	DeleteObject (result1);
	DeleteObject (result2);
	return (_PMathObj)argVal.makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _Constant::IBeta (_PMathObj arg1, _PMathObj arg2) 
{
	if (theValue<=0.0)
	{
		if (theValue < 0.0)
		{
			_String 	errMsg;
			errMsg = _String ("IBeta is defined for x betweeen 0 and 1. Had: ") & theValue;
			ReportWarning 	(errMsg);
		}
		return new _Constant (0.0);
	}

	if (theValue>=1.0)
	{
		if (theValue>1.0)
		{
			_String 	errMsg;
			errMsg = _String ("IBeta is defined for x betweeen 0 and 1. Had: ") & theValue;
			ReportWarning 	(errMsg);
		}
		return new _Constant (1.0);
	}

		
	if ((arg1->ObjectClass()!=NUMBER)||(arg2->ObjectClass()!=NUMBER))
	{
		_String 	errMsg ("IBeta called with a non-scalar argument.");
		WarnError 	(errMsg);
		return 		nil;
	}
	
	_Constant 		*ga = (_Constant*)arg1->Gamma(),
					*gb = (_Constant*)arg2->Gamma();
					
	if (ga&&gb)
	{
		_Constant 	*ac	= (_Constant*)arg1,
					*bc = (_Constant*)arg2;
					
		_Parameter	a = ac->Value(),
					b = bc->Value(),
					x = theValue,
					aa,
					c,
					d,
					del,
					h,
					qab,
					qam,
					qap,
					FPMIN = 1e-100;
					
		bool		swap = false;
					
		long		m,
					m2;
					
		if (x >= (a+1.)/(a+b+2.))
		{
			swap = true;
			c = b;
			b = a;
			a = c;
			x = 1. - x;
		}
		
		qab = a+b;
		qap = a+1.;
		qam = a-1.;
		c   = 1.;
		d   = 1. - qab*x/qap;
		if 	((d<FPMIN)&&(d>-FPMIN))	d = FPMIN;
		d	= 1./d;
		h	= d;
		
		for (m=1;m<100;m++)
		{
			m2 = 2*m;
			aa = m*(b-m)*x / ((qam+m2)*(a+m2));
			d = 1.+aa*d;
			if 	((d<FPMIN)&&(d>-FPMIN))	d = FPMIN;
			c = 1.+aa/c;
			if 	((c<FPMIN)&&(c>-FPMIN))	c = FPMIN;
			d = 1./d;
			h*= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d = 1.+aa*d;
			if 	((d<FPMIN)&&(d>-FPMIN))	d = FPMIN;
			c = 1.+aa/c;
			if 	((c<FPMIN)&&(c>-FPMIN))	c = FPMIN;
			d = 1./d;
			del = d*c;
			h*= del;
			del -= 1.;
			if 	((del<1.e-14)&&(del>-1.e-14))	break;
		}
		
		_Constant 	* res = new _Constant (a+b);
		ac	= (_Constant*)res->Gamma();
		c   = ac->Value()/(ga->Value()*gb->Value()) *
			  exp (a*log(x)+b*log(1-x));
		
		if (swap)
			res->theValue = 1.-c*h/a;
		else
			res->theValue = c*h/a;

		DeleteObject (ac);
		DeleteObject (ga);
		DeleteObject (gb);
		return 	res;
	}
	DeleteObject (ga);
	DeleteObject (gb);
	return nil;
}


//__________________________________________________________________________________
_PMathObj _Constant::IGamma (_PMathObj arg) 
{
	if (arg->ObjectClass()!=NUMBER)
	{
		_String errMsg ("A non-numerical argument passed to IGamma(a,x)");
		WarnError (errMsg);
		return new _Constant (0.0);
	}
	_Parameter x = ((_Constant*)arg)->theValue, sum=0.0;
	if (x>1e25) 
		x=1e25;
	else
		if (x<0)
		{
			_String errMsg ("The domain of x is {x>0} for IGamma (a,x)");
			WarnError (errMsg);
			return new _Constant (0.0);
		}
		else
			if (x==0.0)
				return new _Constant (0.0);
		
	
	if (x<=theValue+1) // use the series representation
	// IGamma (a,x)=exp(-x) x^a \sum_{n=0}^{\infty} \frac{\Gamma((a)}{\Gamma(a+1+n)} x^n
	{
		_Parameter term = 1.0/theValue, den = theValue+1;
		long count = 0;
		while ((fabs(term)>=fabs(sum)*machineEps)&&(count<500))
		{
			sum+=term;
			term*=x/den;
			den += 1.0;
			count++;
		}
	}
	else // use the continue fraction representation
	// IGamma (a,x)=exp(-x) x^a 1/x+/1-a/1+/1/x+/2-a/1+/2/x+...
	{
		_Parameter lastTerm = 0, a0 = 1.0, a1 = x, b0 = 0.0, b1 = 1.0, factor = 1.0, an, ana, anf;
		for (long count = 1; count<500; count++)
		{
			an = count;
			ana = an - theValue;
			a0 = (a1+a0*ana)*factor;
			b0 = (b1+b0*ana)*factor;
			anf = an*factor;
			a1  = x*a0+anf*a1;
			b1  = x*b0+anf*b1;
			if (a1!=0.0)
			{
				factor=1.0/a1;
				sum = b1*factor;
				if (fabs(sum-lastTerm)/sum<machineEps)
					break;
				lastTerm = sum;
			}
		
		}
	}
	_Constant *result = (_Constant*)Gamma();
	result->SetValue(sum*exp(-x+theValue*log(x))/result->theValue);
	if (x>theValue+1)
		result->SetValue (1.0-result->theValue);
	return result;
}

//__________________________________________________________________________________
_PMathObj _Constant::Erf (void) 
{
	_Parameter lV = theValue;
	_Constant  half (.5), sq = (lV*lV);
	_PMathObj  IG = half.IGamma(&sq);
	lV = ((_Constant*)IG)->theValue;
	if (theValue<0)
		lV=-lV;
	((_Constant*)IG)->SetValue(lV);
	return (_PMathObj)IG;
}

//__________________________________________________________________________________
_PMathObj _Constant::ZCDF (void) 
{
	_Parameter lV = theValue;
	
	_Constant  half (.5), 
			   sq (lV*lV/2);
			   
	_PMathObj  IG = half.IGamma(&sq);
	lV = ((_Constant*)IG)->theValue/2;
	
	if (theValue>0)
		((_Constant*)IG)->SetValue(lV+.5);
	else
		((_Constant*)IG)->SetValue(.5-lV);
	return (_PMathObj)IG;
}

//__________________________________________________________________________________
_PMathObj _Constant::Time (void) 
{
	_Constant result;
	if (theValue<1.0)
		result.theValue = ((_Parameter)clock()/CLOCKS_PER_SEC);
	else
	{	
		time_t tt;
		result.theValue = ((_Parameter)time(&tt));
	}
	return 	   (_PMathObj)result.makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _Constant::Less (_PMathObj theObj) 
{
	if (theObj)
		return new _Constant (theValue<((_Constant*)theObj)->theValue);	
	else
		return nil;
}

//__________________________________________________________________________________
_PMathObj _Constant::Greater (_PMathObj theObj) 
{
	if (theObj)
		return new _Constant (theValue>((_Constant*)theObj)->theValue);	
	else
		return nil;
}

//__________________________________________________________________________________
_PMathObj _Constant::GammaDist (_PMathObj alpha, _PMathObj beta) 
{
	_Parameter x = theValue, a = ((_Constant*)alpha)->theValue, 
				b = ((_Constant*)beta)->theValue, gd = exp(a * log(b) -b*x +(a-1)*log(x));
	_Constant * c = (_Constant*)alpha->Gamma();
	gd/=c->theValue;
	c->SetValue(gd);
	return c;
}

//__________________________________________________________________________________
_PMathObj _Constant::CGammaDist (_PMathObj alpha, _PMathObj beta) 
{
	_Parameter     arg = theValue*((_Constant*)beta)->theValue;
	/*if (arg==0)
	{
		_Constant zer (0);
		return    (_PMathObj)zer.makeDynamic();
	}*/
	_Constant newX (arg);
	return alpha->IGamma( &newX);
}

//__________________________________________________________________________________
_PMathObj _Constant::CChi2 (_PMathObj n) 
// chi^2 n d.f. probability up to x
{
	_Constant halfn (((_Constant*)n)->theValue/2), halfx = (theValue/2);
	if ((theValue<0)||(halfn.theValue<0))
	{
		_String warnMsg ("CChi2(x,n) only makes sense for both arguments positive");
		ReportWarning (warnMsg);
		return new _Constant (0.0);
	}
	return halfn.IGamma( &halfx);
}

//__________________________________________________________________________________
_PMathObj _Constant::InvChi2 (_PMathObj n) 
// chi^2 n d.f. probability up to x
{
	if (!chi2)
	{
		_String fla ("IGamma(_n_,_x_)");
		chi2 = new _Formula (fla, nil);
		fla = "_x_^(_n_-1)/Gamma(_n_)/Exp(_x_)";
	    derchi2 = new _Formula (fla,nil);
	}
	_Constant halfn (((_Constant*)n)->theValue*.5);
	if ((theValue<0)||(halfn.theValue<0)||(theValue>1.0))
	{
		_String warnMsg ("InvChi2(x,n) only makes sense for n positive, and x in [0,1]");
		ReportWarning (warnMsg);
		return new _Constant (0.0);
	}
	LocateVar(dummyVariable2)->SetValue (&halfn);
	halfn.SetValue(chi2->Newton(*derchi2,theValue,1e-25,1.e100,LocateVar(dummyVariable1))*2);
	return (_PMathObj)halfn.makeDynamic();	
}

//__________________________________________________________________________________
_PMathObj _Constant::LessEq (_PMathObj theObj) 
{
	if (theObj)
		return new _Constant (theValue<=((_Constant*)theObj)->theValue);	
	else
		return nil;
}

//__________________________________________________________________________________
_PMathObj _Constant::GreaterEq (_PMathObj theObj) 
{
	if (theObj)
		return new _Constant (theValue>=((_Constant*)theObj)->theValue);	
	else
		return nil;
}
//__________________________________________________________________________________
_PMathObj _Constant::AreEqual (_PMathObj theObj) 
{
	if (!theObj) return nil;
	
	_Parameter a = theValue, 
			   b = ((_Constant*)theObj)->theValue;
			   
	if (a==0.0) 
		return new _Constant (b==0.0);
		
	return new _Constant(fabs ((a-b)/a)<tolerance);
}
//__________________________________________________________________________________
_PMathObj _Constant::NotEqual (_PMathObj theObj) 
{
	if (!theObj) return nil;
	_Parameter 	 a = theValue, 
				 b = ((_Constant*)theObj)->theValue;
				 
	if (a==0) 
		return new _Constant (b!=0.0);
		
	return new _Constant(fabs ((a-b)/a)>=tolerance);
}
//__________________________________________________________________________________
_PMathObj _Constant::LAnd (_PMathObj theObj) 
{
	if (!theObj) return nil;
	return new _Constant ((long)(theValue)&&(long)(((_Constant*)theObj)->theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::LOr (_PMathObj theObj) 
{
	if (!theObj) return nil;
	return new _Constant ((long)(theValue)||(long)(((_Constant*)theObj)->theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::LNot () 
{
	return new _Constant (!(long)(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Min (_PMathObj theObj) 
{
	if (!theObj) return nil;
	if (theValue<((_Constant*)theObj)->theValue)
		return (_PMathObj) makeDynamic();
	return 	 (_PMathObj) theObj->makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _Constant::Max (_PMathObj theObj) 
{
	if (!theObj) return nil;
	if (theValue>((_Constant*)theObj)->theValue)
		return (_PMathObj) makeDynamic();
	return 	 (_PMathObj) theObj->makeDynamic();
}

//__________________________________________________________________________________
	
_Variable::_Variable (void)
{
	//hasBeenChanged = false;
	varFormula = nil;
	varFlags = 0;
	theName = nil;
	varValue = nil;
	theIndex = -1;
	SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
}

//__________________________________________________________________________________
void _Variable::Initialize (void)
 {
 	//_Formula::Initialize();
 	_Constant::Initialize();
 	theName = new _String();
 	checkPointer(theName);
 	varValue = nil;
 	theIndex = -1;
	//hasBeenChanged = false;
	varFlags = 0;
	SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
 }
 
 //__________________________________________________________________________________
 void _Variable::Duplicate (BaseRef r)
 {
 	_Variable *v = (_Variable*)r;
	//theFormula.Duplicate (&(v->theFormula));
	if (v->varFormula)
	{
		varFormula = new _Formula();
		varFormula->theFormula.Duplicate (&v->varFormula->theFormula);
	}
	else
		varFormula = nil;
		
	//theStack.theStack.Duplicate (&(v->theStack.theStack));
	//theStack.theStack.Clear();
	theValue = v->theValue;
 	varValue = v->varValue;
 	if (varValue)
 		varValue->nInstances++;
 	theIndex = v->theIndex;
 	theName = v->theName;
 	theName->nInstances++;
 	lowerBound = v->lowerBound;
 	upperBound = v->upperBound;
	//hasBeenChanged = v->hasBeenChanged;
	varFlags = v->varFlags;
 }
 
 //__________________________________________________________________________________
BaseRef _Variable::makeDynamic (void)
{
	_Variable * res = new _Variable;
	if (!res)
	{
		isError(0);
		return nil;
	}
	//memcpy ((char*)res, (char*)this, sizeof (_Variable));
	res->Duplicate(this);
	return res;
}
 
//__________________________________________________________________________________
bool _Variable::CheckFForDependence (long idx, bool opt)
{
	if (varFormula)
		return varFormula->CheckFForDependence (idx, opt);
		
	return false;
}
 
//__________________________________________________________________________________
BaseRef _Variable::toStr(void)
{
	_String res;
	if (varValue&&varValue->IsPrintable())
		return varValue->toStr();
	_PMathObj vv = Compute();
	if (!vv)
		res = "NAN";
	else
	{
		_String* vs = (_String*)vv->toStr();
		res=*vs;
		DeleteObject(vs);
	}
	return res.makeDynamic();
}

//__________________________________________________________________________________
void _Variable::toFileStr(FILE* f)
{
	_String res;
	if (varValue&&varValue->IsPrintable())
	{
		varValue->toFileStr(f);
		return;
	}
	_PMathObj vv = Compute();
	if (!vv)
		fprintf(f,"NAN");
	else
	{
		vv->toFileStr(f);
	}
}
//__________________________________________________________________________________
	
_Variable::_Variable (_String&s, bool isG)
{
	theName = new _String(s);
	checkPointer(theName);
	//hasBeenChanged = false;
	//isGlobal = isG;
	varFlags = isG?HY_VARIABLE_GLOBAL:0;
	varValue = nil;
	varFormula = nil;
	SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
	InsertVar (this);
}

//__________________________________________________________________________________
	
_Variable::_Variable (_String&s, _String&f, bool isG)//:  _Formula (f)
{
	//hasBeenChanged = false;
	//isGlobal = isG;
	theName = new _String(s);
	checkPointer(theName);
	varFlags = isG?HY_VARIABLE_GLOBAL:0;
	varValue = nil;
	SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
	InsertVar (this);
	varFormula = new _Formula (f);
	if (varFormula->IsAConstant())
	{
		_PMathObj theP =varFormula->Compute();
		if (theP)
		{
			SetValue (theP);
			delete (varFormula);
			varFormula = nil;
			//theFormula.Clear();
			//theStack.theStack.Clear();
		}
		else
			return;
	}
}



//__________________________________________________________________________________
	
_Variable::~_Variable (void)
{
	nInstances++;
	if (varValue)
	{
		DeleteObject (varValue);
	}
	if (theName)
	{
		DeleteObject (theName);
	}
	if (varFormula)
		delete (varFormula);
}
	
//__________________________________________________________________________________
bool	_Variable::IsVariable (void) {return true;}

//__________________________________________________________________________________

_PMathObj  _Variable::Compute (void) // compute or return the value 
{
	if (varFormula == nil) // no formula, just return the value
	{
		if (varValue)
			return varValue->Compute();

		varValue =  new _Constant(theValue);
	}
	else
	{
		if (useGlobalUpdateFlag)
		{			
			if ((varFlags & HY_DEP_V_COMPUTED) && varValue)
				return varValue;
			else
				if (varFormula->HasChanged()||!varValue)
				{
					DeleteObject (varValue);
					varValue = (_PMathObj)varFormula->Compute()->makeDynamic();
				}

			varFlags |= HY_DEP_V_COMPUTED;
		}
		else
			if (varFormula->HasChanged()||!varValue)
			{
				DeleteObject (varValue);
				varValue = (_PMathObj)varFormula->Compute()->makeDynamic();
			}

	}

	return varValue;
}

//__________________________________________________________________________________

void  _Variable::CompileListOfDependents (_SimpleList& rec) 
{
	_SimpleList tcache;
	long		iv,
				i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());
				
	for (; i >= 0; i = variableNames.Traverser (tcache,iv))
	{
		_Variable* thisVar = FetchVar (i);
		if (!thisVar->IsIndependent())
		{
			if (thisVar->CheckFForDependence (theIndex))
			{
				long f = thisVar->GetAVariable();
				if (rec.Find(f)<0)
					rec<<f;
			}
		}
	}
}

//__________________________________________________________________________________
void  _Variable::SetValue (_PMathObj theP, bool dup) // set the value of the var
{
	//hasBeenChanged = true;
	varFlags |= HY_VARIABLE_CHANGED;
	if (theP->ObjectClass()==NUMBER)
	{
		if (varFormula)
		{
		 
			// also update the fact that this variable is no longer dependent in all declared
			// variable containers which contain references to this variable
			long i;
			for (i = 0; i<variablePtrs.lLength; i++)
			{
				if (freeSlots.Find(i)>=0) continue;
				_Variable* theV = (_Variable*)variablePtrs(i);
				if (theV->IsContainer())
				{
					_VariableContainer* theVC = (_VariableContainer*)theV;
					if (!theVC->RemoveDependance (theIndex))
					{
					   ReportWarning ((_String("Can't make variable ")&*GetName()&" independent in the context of "&*theVC->GetName()&" because its template variable is not independent."));
					   continue;
					}   
				}
			}
			for (i = 0; i<likeFuncList.lLength; i++)
				if (((_String*)likeFuncNamesList(i))->sLength)
					((_LikelihoodFunction*)likeFuncList(i))->UpdateDependent(theIndex);
				
			//_Formula::Clear();
			delete (varFormula);
			varFormula = nil;
		}
		if (varValue)
		{	
			DeleteObject (varValue);
			varValue=nil;
		}
		
		theValue = theP->Value();
		
		if (!dup)
			DeleteObject (theP);
			
		if ((theValue<lowerBound)||(theValue>upperBound))
		{
			/*_String *t = (_String*)theP->toStr(),s;
			s = _String("Value:")&*t&" is out of bounds for the variable "&*theName;
			ReportWarning (s);
			DeleteObject(t);*/
			if (theValue<=lowerBound+1e-50)
				theValue = lowerBound;
			else
				theValue = upperBound;
		}
		
	}
	else
	{
		if (varFormula)
		{
			delete (varFormula);
			varFormula = nil;
			//theFormula.Clear();
		}
		if (varValue)
		{	
			DeleteObject (varValue);
			varValue=nil;
		}
		if (theP->ObjectClass()==TREE)
		{
			variablePtrs.lData[theIndex] = (long)(((_TheTree*)theP)->makeDynamicCopy(GetName()));
			DeleteObject(this);
		}
		else
			if (dup)
				varValue = (_PMathObj)theP->makeDynamic();
			else
				varValue = theP;			
	}
}

//__________________________________________________________________________________
void  _Variable::SetNumericValue (_Parameter v) // set the value of the var to a number
{
	//hasBeenChanged = true;
	varFlags |= HY_VARIABLE_CHANGED;
	theValue = v;
			
	if (theValue<lowerBound || theValue>upperBound)
	{
		if (theValue<=lowerBound+1e-50)
			theValue = lowerBound;
		else
			theValue = upperBound;
	}
}
		
//__________________________________________________________________________________

void  _Variable::CheckAndSet (_Parameter c, bool oob) // set the value of the var
{
	//hasBeenChanged = true;
	varFlags |= HY_VARIABLE_CHANGED;
	_Parameter l = lowerBound+1.0e-30,
			   u = upperBound-1.0e-30;
	if (c<l || c>u )
	{
		if (oob)
			return;
			
		if (c<l)
		{
			theValue = l;
		}
		else
			theValue = u;
	}
	else
	{
		theValue = c;
	}

	if (varValue)
		DeleteObject (varValue);

	varValue =  new _Constant(theValue);
}

//__________________________________________________________________________________
void	_Variable::ClearConstraints (void)
{
	if (IsCategory ())
	{
		_Variable newVar (*GetName(), IsGlobal());
		newVar.SetValue ((_PMathObj)Compute()->makeDynamic(),false);
		ReplaceVar ( &newVar);
		/*_Matrix * modelMatrix = (_Matrix*)LocateVar(modelMatrixIndices.lData[1])->GetValue();
		for (long k=0; k<4; k++)
			for (long k2 = 0; k2<4; k2++)
				if (k!=k2)
				{
					StringToConsole (*(_String*)modelMatrix->GetFormula(k,k2)->toStr());
					BufferToConsole ("\n");
				}
		*/
	}
	else
	{
		if (!IsIndependent()) 
			SetValue ((_PMathObj)Compute()->makeDynamic(),false);
		SetBounds (DEFAULTLOWERBOUND,DEFAULTUPPERBOUND);
	}
}

//__________________________________________________________________________________
bool _Variable::IsConstant (void)
{ 
	if (varFormula && varFormula->theFormula.lLength)
		return varFormula->IsConstant();
	
	if (varValue && varValue->ObjectClass () != NUMBER)
		return varValue->IsConstant();
	
	return false;
}

//__________________________________________________________________________________

void  _Variable::SetFormula (_Formula& theF) // set the value of the var to a formula
{
	bool changeMe    = false,
		 isAConstant = theF.IsAConstant();
	
	_Formula* myF = &theF;
	
	if (isAConstant)
	{
		_PMathObj theP = theF.Compute();
		if (theP)
		{
			myF = new _Formula ((_PMathObj)theP->makeDynamic(),false);
			checkPointer (myF);
		}
		else
			return;
	}
		
	_SimpleList vars;
	{
		_AVLList vA (&vars);
		theF.ScanFForVariables (vA,true);
		vA.ReorderList();
	}
	
	if (vars.BinaryFind(theIndex)>=0)
	{
		_String * sf = (_String*)theF.toStr();
		WarnError ((_String("Can't set variable ")&*GetName()&" to "&*sf&" because it would create a circular dependance."));
		DeleteObject(sf);
		if (&theF!=myF)
			delete myF;
		return;
	}
	
	
	if (varFlags & HY_VARIABLE_CHANGED)
		varFlags -= HY_VARIABLE_CHANGED;


	if (varFormula)
	{
		delete (varFormula);
		varFormula = nil;
	}
	else
	  	changeMe = true;
	  
	if (varValue)
	{	
		DeleteObject (varValue);
		varValue=nil;
	}

	//_Formula::Duplicate ((BaseRef)myF);
	varFormula = new _Formula;
	varFormula->Duplicate ((BaseRef)myF);
	
	// mod 20060125 added a call to simplify constants
	varFormula->SimplifyConstants ();
	
	// also update the fact that this variable is no longer independent in all declared
	// variable containers which contain references to this variable
	if (changeMe)
		if (deferSetFormula)
		{
			*deferSetFormula << theIndex;
			deferIsConstant  << isAConstant;
		}
		else
		{
			long i;
			_SimpleList tcache;
			long		iv;
			
			i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());
						
			for (; i >= 0; i = variableNames.Traverser (tcache,iv))
			{
				_Variable* theV = FetchVar(i);
				if (theV->IsContainer())
				{
					_VariableContainer* theVC = (_VariableContainer*)theV;
					if (theVC->SetDependance(theIndex) == -2)
					{
				       ReportWarning ((_String("Can't make variable ")&*GetName()&" dependent in the context of "&*theVC->GetName()&" because its template variable is bound by another relation in the global context."));
				       continue;
				    } 
				}
			}
			{
			for (long i = 0; i<likeFuncList.lLength; i++)
				if (((_String*)likeFuncNamesList(i))->sLength)
					((_LikelihoodFunction*)likeFuncList(i))->UpdateIndependent(theIndex,isAConstant);
        	}
		}
	
	if (&theF!=myF)
		delete myF;
}

//__________________________________________________________________________________
void  _Variable::PreMarkChanged (void)
{
	if (varFormula)	
	{
		varFlags &= HY_DEP_V_INSPECTED_CLR;
		
		if (HasChanged(false))
			varFlags |= HY_DEP_V_MODIFIED;
		if (HasChanged(true))
			varFlags |= HY_DEP_V_MODIFIED_CATS;

		varFlags |= HY_DEP_V_INSPECTED;
	}	
}

//__________________________________________________________________________________
void  _Variable::PostMarkChanged (void)
{
	varFlags &= HY_DEP_CLEAR_MASK;
}


//__________________________________________________________________________________
bool  _Variable::HasChanged (bool ignoreCats) // does this variable need recomputing
{
	if (varFormula)	
	{
		if (useGlobalUpdateFlag && (varFlags&HY_DEP_V_COMPUTED))
			return false;
		
		if (varFlags&HY_DEP_V_INSPECTED)
			return ignoreCats?(varFlags&HY_DEP_V_MODIFIED_CATS):(varFlags&HY_DEP_V_MODIFIED);
		
		return	varFormula->HasChanged(ignoreCats);

	}
	else
	{
		if (varValue&&(varValue->IsVariable()))
			return varValue->HasChanged();
		if (ignoreCats && IsCategory())
			return false;
		return varFlags & HY_VARIABLE_CHANGED;
	}
	
}

//__________________________________________________________________________________
				
void _Variable::MarkDone (void)
{
	if (!varFormula && (varFlags & HY_VARIABLE_CHANGED) && !(varValue && varValue->IsVariable()))
		varFlags -= HY_VARIABLE_CHANGED;
}

//__________________________________________________________________________________

_VariableContainer::_VariableContainer (void)
{
	theParent = nil;
	theModel = -1;
	iVariables = nil;
	dVariables = nil;
	gVariables = nil;
}

//__________________________________________________________________________________

void	_VariableContainer::Duplicate (BaseRef theO)
{
	_Variable::Duplicate (theO);
	_VariableContainer *theVC = (_VariableContainer*)theO;
	theParent= theVC->theParent;
	theModel = theVC->theModel;
	if (theVC->iVariables)
	{
		if (iVariables)
			iVariables->Clear();
		else
			checkPointer(iVariables = new _SimpleList);
		iVariables->Duplicate (theVC->iVariables);
	}
	else
	{
		if (iVariables)
		{
			delete (iVariables);
			iVariables = nil;
		}
	}
	if (theVC->dVariables)
	{
		if (dVariables)
			dVariables->Clear();
		else
			checkPointer(dVariables = new _SimpleList);
		dVariables->Duplicate (theVC->dVariables);
	}
	else
	{
		if (dVariables)
		{
			delete (dVariables);
			dVariables = nil;
		}
	}
	if (theVC->gVariables)
	{
		if (gVariables)
			gVariables->Clear();
		else
			checkPointer (gVariables = new _SimpleList);
		gVariables->Duplicate (theVC->gVariables);
	}
	else
	{
		if (gVariables)
		{
			delete (gVariables);
			gVariables = nil;
		}
	}
}


//__________________________________________________________________________________

void	_VariableContainer::TrimMemory ()
{
	if (iVariables)
		iVariables->TrimMemory();
	if (dVariables)
		dVariables->TrimMemory();
	if (gVariables)
		gVariables->TrimMemory();
}

//__________________________________________________________________________________

BaseRef _VariableContainer::makeDynamic (void)
{
	_VariableContainer * res = new _VariableContainer;
	checkPointer(res);
	memcpy ((char*)res, (char*)this, sizeof (_VariableContainer)); // ???
	res->Duplicate(this);
	return res;
}

//__________________________________________________________________________________

BaseRef _VariableContainer::toStr (void)
{
	_String * res = new _String (128L,true);
	
	checkPointer (res);
	
	(*res) << "Container Class:";
	(*res) << theName;
	(*res) << ":{ Independent Variables:";
	
	if (iVariables)
		for (long i = 0; i<iVariables->lLength; i+=2)
		{
			_String* s = (_String*)variablePtrs(iVariables->lData[i])->toStr();
			(*res) << s;
			#ifndef USE_POINTER_VC
			if (i<independentVars.lLength-1) 
			#else
			if (i<iVariables->lLength-2) 		
			#endif
				(*res) << ',';
			DeleteObject(s);
		}
	
	(*res) << "; Dependent Variables:";
	
	if (dVariables)
		for (long i2 = 0; i2<dVariables->lLength; i2+=2)
		{
			_String* s = (_String*)variablePtrs(dVariables->lData[i2])->toStr();
			(*res) << s;
			#ifndef USE_POINTER_VC
			if (i2<independentVars.lLength-1) 
			#else
			if (i2<dVariables->lLength-2) 		
			#endif
				(*res) << ',';
			DeleteObject(s);
		}
	
	(*res) << '}';
	res->Finalize();
	return res;
}

//__________________________________________________________________________________

_VariableContainer::_VariableContainer (_String theName, _String theTmplt, _VariableContainer* theP)
{
	iVariables = nil;
	dVariables = nil;
	gVariables = nil;
	InitializeVarCont (theName, theTmplt,theP);
}

//__________________________________________________________________________________

bool _VariableContainer::HasExplicitFormModel (void)
{
	if (theModel == -1)
		return false;
	return (modelTypeList.lData[theModel]);
}


//__________________________________________________________________________________

_Matrix* _VariableContainer::GetModelMatrix (void)
{
	if (theModel == -1)
		return nil;
	if (modelTypeList.lData[theModel])
		return (_Matrix*) ((_Formula*)modelMatrixIndices.lData[theModel])->Compute();
	return (_Matrix*) (LocateVar(modelMatrixIndices.lData[theModel])->GetValue());
}

//__________________________________________________________________________________

long _VariableContainer::GetModelDimension (void)
{
	long matrixDim = 0;
	if (theModel >= 0)
	{
		matrixDim = modelTypeList.lData[theModel];
		if (matrixDim == 0)
			return GetModelMatrix()->GetHDim();
	}
	return matrixDim;
}

//__________________________________________________________________________________

_Matrix* _VariableContainer::GetFreqMatrix (void)
{
	if (theModel>=0)
	{
		long freqID = modelFrequenciesIndices.lData[theModel];
		if (freqID>=0)
			return (_Matrix*) (LocateVar(freqID)->GetValue());
		else
			return (_Matrix*) (LocateVar(-freqID-1)->GetValue());
	}
	return nil;
}

//__________________________________________________________________________________
void	_VariableContainer::InitializeVarCont (_String& aName, _String& theTmplt, _VariableContainer* theP, _AVLListXL* varCache)
{
	_String fullName (aName), 
	varName;

	theModel  = FindModelName(theTmplt);
	theParent = theP;
	
	
	if (aName.sLength)
	{
		long f = aName.Find('.');
		
		while (theP)
		{
			if (f!=-1)
				f = aName.Find('.',f+1,-1);
			else
				break;
				
			theP = theP->theParent;
		}
		
		if (theP) 			
			fullName = (*(theP->theName))&'.'&fullName;

		theName = (_String*)(fullName.makeDynamic());
	    InsertVar (this);
    }
    else
    	fullName = *theName;
    	
   
	if (theModel!= HY_NO_MODEL) // build the matrix variables
	{
		_SimpleList 	  mVars;
		{
			long cachedID = -1;
			bool doScan   = !varCache || (cachedID = varCache->Find ((BaseRef) theModel)) < 0 ;
			
			if (doScan)
			{
			
				_AVLList				ma (&mVars);
				ScanModelForVariables   (GetModelIndex(), ma,true,theModel,false);
				
				long freqID 	= modelFrequenciesIndices.lData[theModel];
				if (freqID>=0)
					((_Matrix*) (LocateVar(freqID)->GetValue()))->ScanForVariables2(ma,true,-1,false);
					
				ma.ReorderList();
				
				if (varCache)
					varCache->Insert ((BaseRef)theModel, (long)mVars.makeDynamic(),false);
			}
			else
				mVars.Duplicate (varCache->GetXtra (cachedID));

		}
		
		for (long i=0; i<mVars.lLength; i++)
		{
			_Variable * aVar = (_Variable*)variablePtrs (mVars.lData[i]);
			if (aVar->IsGlobal()) 
			{
				//if (curVar->IsIndependent())
				{
					if (!gVariables)
						checkPointer (gVariables = new _SimpleList);
					(*gVariables) << aVar->GetAVariable();
				}
				continue;
			}
			
			long f = aVar->theName->FindBackwards('.',0,-1);
			if (f>=0)
				varName = fullName&'.'& aVar->theName->Cut(f+1,-1);
			else
				varName = fullName&'.'& *aVar->theName;
			
			 
			f = LocateVarByName (varName);
			if (f<0)
			{
				_Variable v (varName);
				f = v.theIndex;
			}
			else
				f = variableNames.GetXtra (f);
			
			_Variable * spawnedVar = FetchVar (f);
			spawnedVar->SetBounds (aVar->GetLowerBound(), aVar->GetUpperBound());

			if (aVar->IsIndependent())
			{
				if (!iVariables)
					checkPointer (iVariables = new _SimpleList);
				(*iVariables) << f;
				(*iVariables) << mVars.lData[i];
			}
			else
			{
				if (!dVariables)
					checkPointer (dVariables = new _SimpleList);
				(*dVariables) << f;
				(*dVariables) << mVars.lData[i];
			}
		}
	}
	SortVars();
}

//__________________________________________________________________________________

void _VariableContainer::ScanAndAttachVariables (void)
{
	_Variable* curVar;
	_SimpleList travcache;

	long f = variableNames.Find (theName,travcache);
	if (f<0) 
		return;
	
	_String theNameAndADot (*theName);
	theNameAndADot = theNameAndADot&'.';
	
	
			
	for (f = variableNames.Next (f, travcache);f>=0;f = variableNames.Next (f, travcache))
	{
		curVar = FetchVar (f);
		
		if (curVar->theName->startswith(theNameAndADot))
		{
			if (!curVar->IsContainer()) 
			{
				long   vix = variableNames.GetXtra (f);
				
				if (curVar->IsIndependent())
				{
					if ( ((!iVariables)||iVariables->FindStepping(vix,2)==-1) && ((!dVariables)||dVariables->FindStepping(vix,2)==-1)) 
					{
						if (!iVariables)
							checkPointer (iVariables = new _SimpleList);
						(*iVariables)<<vix;
						(*iVariables)<<-1;
					}						
				}
				else
				{
					if ( ((!iVariables)||iVariables->FindStepping(vix,2)==-1) && ((!dVariables)||dVariables->FindStepping(vix,2)==-1)) 
					{
						if (!dVariables)
							checkPointer (dVariables = new _SimpleList);
						(*dVariables)<<vix;
						(*dVariables)<<-1;
					}						
				}
			}
		}
		else	
			break;
	}
	


}
//__________________________________________________________________________________
				
_VariableContainer::~_VariableContainer(void)
{
	if (iVariables) delete iVariables;
	if (dVariables) delete dVariables;
	if (gVariables) delete gVariables;
}
	
//__________________________________________________________________________________
				
bool _VariableContainer::HasChanged (void)
{
	long i;
	if (iVariables)
		for (i = 0; i<iVariables->lLength; i+=2)
			if (LocateVar (iVariables->lData[i])->HasChanged()) return true;
	
	if (gVariables)
		for (i = 0; i<gVariables->lLength; i++)
			if (LocateVar (gVariables->lData[i])->HasChanged()) return true;
	
	if (dVariables)
		for (i = 0; i<dVariables->lLength; i+=2)
			if (LocateVar (dVariables->lData[i])->HasChanged()) return true;
	
	return false;
}

//__________________________________________________________________________________
				
_Variable* _VariableContainer::GetIthIndependent (long index)
{
	if (iVariables && (index*=2)<iVariables->lLength)
		return LocateVar (iVariables->lData[index]);
	else 	
		return nil;
}

//__________________________________________________________________________________
				
_Variable* _VariableContainer::GetIthDependent (long index)
{
	if (dVariables && (index*=2)<dVariables->lLength)
		return LocateVar (dVariables->lData[index]);
	else 	
		return nil;
}

//__________________________________________________________________________________
				
_Variable* _VariableContainer::GetIthParameter (long index)
{		
	if (iVariables)
	{
		if ( (index*=2) <iVariables->lLength)
			return LocateVar (iVariables->lData[index]);
		else 	
		{
			if (dVariables)
			{
				index-=iVariables->lLength;
				if (index<dVariables->lLength)
					return LocateVar (dVariables->lData[index]);
			}
		}
	}
	else
	{
		if (dVariables && (index*=2) <dVariables->lLength)
			return LocateVar (dVariables->lData[index]);
	}
	return nil;
}

//__________________________________________________________________________________
				
bool _VariableContainer::NeedToExponentiate (bool ignoreCats)
{
	if (HY_VC_NO_CHECK&varFlags)
		return false;
	
	if (iVariables)
		for (long i = 0; i<iVariables->lLength && iVariables->lData[i+1] >= 0; i+=2)
		{
			if (LocateVar (iVariables->lData[i])->HasChanged(ignoreCats)) 
				return true;
		}
	
	if (gVariables)
		for (long i = 0; i<gVariables->lLength; i++)
			if (LocateVar (gVariables->lData[i])->HasChanged(ignoreCats)) 
				return true;
	if (dVariables)
		for (long i = 0; i<dVariables->lLength && dVariables->lData[i+1] >= 0; i+=2)
			if (LocateVar (dVariables->lData[i])->HasChanged(ignoreCats)) 
				return true;
	
	return false;
}

//__________________________________________________________________________________
void	  _VariableContainer::SortVars(void)
{
	// sort independents 1st
	// use dumb bubble sort
	bool 		done = false;
	long 		t, 
				index;
				
	_String		*s1, 
				*s2;
				
	if (iVariables && iVariables->lLength>2)
	{
		while (!done)
		{
			done = true;
			s1 = LocateVar(iVariables->lData[0])->GetName();
			for (index = 2; index<iVariables->lLength; index+=2)
			{
				s2 = LocateVar(iVariables->lData[index])->GetName();
				if (s2->Less (s1))
				{
					done = false;
					t = iVariables->lData[index];
					iVariables->lData    [index] = iVariables->lData[index-2];
					iVariables->lData    [index-2] = t;
					
					t = iVariables->lData[index+1];
					iVariables->lData    [index+1] = iVariables->lData[index-1];
					iVariables->lData    [index-1] = t;
				}
				
			}
		}
	}
	if (dVariables && dVariables->lLength>2)
	{
		done = false;
		while (!done)
		{
			done = true;
			s1 = LocateVar(dVariables->lData[0])->GetName();
			for (index = 2; index<dVariables->lLength; index+=2)
			{
				s2 = LocateVar(dVariables->lData[index])->GetName();
				if (s2->Less (s1))
				{
					done = false;
					t = dVariables->lData[index];
					dVariables->lData    [index] = dVariables->lData[index-2];
					dVariables->lData    [index-2] = t;
					
					t = dVariables->lData[index+1];
					dVariables->lData    [index+1] = dVariables->lData[index-1];
					dVariables->lData    [index-1] = t;
				}
				
			}
		}
	}
}
//__________________________________________________________________________________
bool	  _VariableContainer::RemoveDependance (long varIndex)
{
	if (dVariables)
	{
		long f = dVariables->FindStepping(varIndex,2);
				 
		if (f!=-1)
		{
			
				//printf ("Moving dep->ind for %s from %s\n", LocateVar (varIndex)->GetName()->sData,
				//		GetName()->sData);
			  
			   
			   /*if (dVariables->lData[f+1]>=0)
			   {
				  _Variable* checkVar = LocateVar(dVariables->lData[f+1]);
				   printf ("Local variable %s\n", checkVar->GetName()->sData);
				  //if (!checkVar->IsIndependent())
					//  return false;
			   }*/
			
		       _String* thisName = LocateVar (dVariables->lData[f])->GetName();
		       
		       long insPos = 0;
			
		       if (!iVariables)
		       		checkPointer (iVariables = new _SimpleList);
		       
		       while (insPos<iVariables->lLength && (thisName->Greater (LocateVar (iVariables->lData[insPos])->GetName())))
		       		insPos+=2;
		       		
			   
			   iVariables->InsertElement ((BaseRef)varIndex, insPos, false, false);
			   iVariables->InsertElement ((BaseRef)dVariables->lData[f+1], insPos+1, false, false);
			   
			   if (dVariables->lLength>2)
			   {
			   		dVariables->Delete(f);
			   		dVariables->Delete(f);
			   		dVariables->TrimMemory();
			   }
			   else
			   {
				   delete dVariables;
				   dVariables = nil;
			   }
		}
	}
	return true;
}

//__________________________________________________________________________________
long	  _VariableContainer::CheckAndAddUserExpression (_String& pName, long startWith)
{
	_String tryName, tryName2;
	tryName = (*theName)&'.'&pName;
	tryName2 = tryName;
	long 	k = startWith>2?startWith:2;
	if (startWith>=2)
		tryName2 = tryName&startWith;
		
	while (LocateVarByName(tryName2)>=0)
	{
		tryName2 = tryName&k;
		k++;
	}
	
	if (startWith<0)
		return k>2?k-1:0;
		
	if (startWith<2)
	{
		if (k>2)
			pName = pName&_String(k-1);
	}
	else
	{
		if (k>startWith)
			pName = pName & _String (k-1);
		else
			pName = pName & _String (startWith);
	}
		
	_Variable newVar (tryName2);
	k =  newVar.GetAVariable();
	
	if (!dVariables)
		checkPointer (dVariables = new _SimpleList);
	(*dVariables) << k;
	(*dVariables) << -1;
	return k;
}

//__________________________________________________________________________________
void	  _VariableContainer::CopyMatrixParameters (_VariableContainer* source)
{
	if (iVariables && source->iVariables)
	for (long i=0; i<iVariables->lLength && i< source->iVariables->lLength;i+=2)
		LocateVar (iVariables->lData[i])->SetValue(LocateVar (source->iVariables->lData[i])->Compute());

	_PMathObj srcVal = source->Compute();
	SetValue (srcVal);
}

//__________________________________________________________________________________
void	  _VariableContainer::KillUserExpression (long varID)
{
	if (dVariables)
	{
		long f = dVariables->FindStepping(varID,2);
		if (f>=0)
		{
			DeleteVariable (*LocateVar(varID)->GetName(),true);
			if (dVariables->lLength > 2)
			{
				dVariables->Delete (f);
				dVariables->Delete (f);
				dVariables->TrimMemory ();
			}
			else
			{
				delete dVariables;
				dVariables = nil;
			}
		}
	}
}

//__________________________________________________________________________________
long	  _VariableContainer::SetDependance (long varIndex)
{
	if (iVariables)
	{
		long f;

		if (varIndex>=0)
		{
			f = iVariables->FindStepping(varIndex,2);
			if (f<0)
				return -1;	
		}
		else
		{
			f = -varIndex-1;
			varIndex = iVariables->lData[f];
		}
	
			
		//printf ("Moving ind->dep for %s from %s\n", LocateVar (varIndex)->GetName()->sData,
		//		GetName()->sData);
		
		if (iVariables->lData[f+1]>=0)
	    {
			//printf ("Local variable %s\n", LocateVar (iVariables->lData[f+1])->GetName()->sData);
			if (!LocateVar(iVariables->lData[f+1])->IsIndependent())
		       return -2;
        }

        _String* thisName = LocateVar (iVariables->lData[f])->GetName();
        
        long	insPos = 0;
        
        if (!dVariables)
        	checkPointer (dVariables = new _SimpleList);
       	
       	while (insPos<dVariables->lLength)
		{
			_Variable *dVar = LocateVar (dVariables->lData[insPos]);
			if (!dVar)
			{
				FlagError ("Internal error in SetDependance()");
				return -1;
			}
			if (!thisName->Greater (dVar->GetName()))
				break;
       		insPos+=2;
		}
		
		dVariables->InsertElement ((BaseRef)varIndex, insPos, false, false);
		dVariables->InsertElement ((BaseRef)iVariables->lData[f+1], insPos+1, false, false);
		
		if (iVariables->lLength > 2)
		{
			iVariables->Delete(f);
			iVariables->Delete(f);
			iVariables->TrimMemory();
		}
		else
		{
			delete iVariables;
			iVariables = nil;
		}
		
		return varIndex;
	}
	return -1;
}

//__________________________________________________________________________________
bool	  _VariableContainer::SetMDependance (_SimpleList& mDep)
{
	if (iVariables)
		if (mDep.lLength*2 > iVariables->lLength)
			for (long k=iVariables->lLength-2; k>=0; k-=2)
			{
				long f = mDep.BinaryFind (iVariables->lData[k]);
				if (f>=0)
					SetDependance (-k-1);
			}
		else
			for (long k=0; iVariables && k<mDep.lLength; k++)
				SetDependance (mDep.lData[k]);
		
	return true;
}


//__________________________________________________________________________________
void	  _VariableContainer::Clear(void)
{
	theModel = -1;
	if (iVariables)
	{
		delete iVariables;
		iVariables = nil;
	}
	if (dVariables)
	{
		delete dVariables;
		dVariables = nil;
	}
	if (gVariables)
	{
		delete gVariables;
		gVariables = nil;
	}
}

//__________________________________________________________________________________
long	  _VariableContainer::CountAll(void)
{
	return (iVariables?iVariables->lLength/2:0)+(dVariables?dVariables->lLength/2:0);
}

//__________________________________________________________________________________
long	  _VariableContainer::CountIndependents(void)
{
	return iVariables?iVariables->lLength/2:0;
}

//__________________________________________________________________________________
bool	  _VariableContainer::HasLocals  (void) 
{
	return (iVariables && iVariables->lLength)||(dVariables&&dVariables->lLength);
}

//__________________________________________________________________________________
bool	  _VariableContainer::IsModelVar  (long i) 
{
	return dVariables->lData[2*i+1]>=0;
}

//__________________________________________________________________________________

_String*	_VariableContainer::GetSaveableListOfUserParameters (void)
{
	_String * result = new _String (64L, true);
	checkPointer (result);
	
	if (dVariables)
		for (long i=0; i<dVariables->lLength; i+=2)
			if (dVariables->lData[i+1]<0)
			{
				_Variable * userParm  = (_Variable*) LocateVar (dVariables->lData[i]);
				_String	  * varString = (_String*)userParm->GetFormulaString();
				*result << userParm->GetName();
				*result << ':';
				*result << '=';
				*result << varString;
				DeleteObject (varString);
				*result << ';';
				*result << '\n';
			}	
	
	result->Finalize();
	return result;
}

//__________________________________________________________________________________
void	  _VariableContainer::ClearConstraints(void)
{
	while (dVariables)
		LocateVar(dVariables->lData[0])->ClearConstraints();
}

//__________________________________________________________________________________

void  _VariableContainer::CompileListOfDependents (_SimpleList& rec) 
{
	if (iVariables)
		for (long i=0; i<iVariables->lLength; i+=2)
			LocateVar(iVariables->lData[i])->CompileListOfDependents (rec);

	if (gVariables)
		for (long i=0; i<gVariables->lLength; i++)
			LocateVar(gVariables->lData[i])->CompileListOfDependents (rec);

	if (dVariables)
	{
		for (long i=0; i<dVariables->lLength; i+=2)
			LocateVar(dVariables->lData[i])->CompileListOfDependents (rec);

		{
		for (long i=0; i<dVariables->lLength; i+=2)
		{
			long f = rec.Find (dVariables->lData[i]);
			if (f>=0)
				rec.Delete (f);
		}
        }
	}
}


//__________________________________________________________________________________
				
void _VariableContainer::MarkDone (void)
{
	if (iVariables)
		for (long i = 0; i<iVariables->lLength && iVariables->lData[i+1] >= 0; i+=2)
			LocateVar (iVariables->lData[i])->MarkDone();
	if (gVariables)
		for (long i = 0; i<gVariables->lLength; i++)
			LocateVar (gVariables->lData[i])->MarkDone();	
}

//__________________________________________________________________________________
				
void _VariableContainer::MatchParametersToList (_List& suffixes, bool doAll, bool indOnly)
{
	if (doAll)
	{
		for (long i=suffixes.lLength-1; i>=0; i--)
		{
			long j;
			if (!indOnly)
			{
				if (dVariables)
				{
					for (j=0; j<dVariables->lLength; j+=2)
						if (LocateVar(dVariables->lData[j])->GetName()->endswith (*(_String*)suffixes.lData[i]))
							break;

					if (j<dVariables->lLength)
						continue;
				}
			}
			if (iVariables)
			{
				for (j=0; j<iVariables->lLength; j+=2)
				{
					if (LocateVar(iVariables->lData[j])->GetName()->endswith (*(_String*)suffixes.lData[i]))
						break;
				}
				if (j==iVariables->lLength)
					suffixes.Delete (i);
			}
			else
				suffixes.Delete (i);
		}
	}
	else
	{
		for (long i=suffixes.lLength-1; i>=0; i--)
		{
			long j;
			if (dVariables)
			{
				for (j=0; j<dVariables->lLength; j+=2)
				{
					if (dVariables->lData[j+1]<0)
					{
						if (LocateVar(dVariables->lData[j])->GetName()->endswith (*(_String*)suffixes.lData[i]))
							break;
					}
				}
				if (j==dVariables->lLength)
					suffixes.Delete (i);
			}
			else
				suffixes.Delete(i);
		}
	}
}

//__________________________________________________________________________________
				
bool _VariableContainer::IsConstant (void)
{
	if (iVariables)
		return false;
	
	if (dVariables)
		for (long i = 0; i<dVariables->lLength; i+=2)
			if (!LocateVar(dVariables->lData[i])->IsConstant())
				return false;
			
	if (gVariables)
		for (long i = 0; i<gVariables->lLength; i++)
			if (!LocateVar(gVariables->lData[i])->IsConstant())
				return false;
			
	return true;
}

//__________________________________________________________________________________
				
void _VariableContainer::ScanForVariables (_AVLList& l,_AVLList& l2)
{
	if (iVariables)
		for (long i = 0; i<iVariables->lLength; i+=2)
			l.Insert ((BaseRef)iVariables->lData[i]);
	if (dVariables)
		for (long i = 0; i<dVariables->lLength; i+=2)
		{	
			l2.Insert ((BaseRef)dVariables->lData[i]);
			_SimpleList temp;
			{
				_AVLList  ta (&temp);
				LocateVar (dVariables->lData[i])->ScanForVariables(ta, true);
				ta.ReorderList();
			}
			// see if any of them are global
			for (long j=0; j<temp.lLength;j++)
			{
				long p = temp.lData[j];
				_Variable * v = LocateVar(p);
				if (!v->IsGlobal() && v->IsIndependent())
					l.Insert ((BaseRef)p);
			}
		}
}

//__________________________________________________________________________________
				
void _VariableContainer::ScanForDVariables (_AVLList& l,_AVLList&)
{
	if (dVariables)
		for (long i = 0; i<dVariables->lLength; i+=2)
			l.Insert ((BaseRef)dVariables->lData[i]);	
}

//__________________________________________________________________________________
				
void _VariableContainer::GetListOfModelParameters (_List& rec)
{
	if (iVariables)
		for (long i = 1; i<iVariables->lLength; i+=2)
		{
			long p = iVariables->lData[i];
			if (p>=0)
				rec << LocateVar(p)->GetName();
		}	
}

//__________________________________________________________________________________
				
void _VariableContainer::ScanForGVariables (_AVLList& l,_AVLList& l2)
{
	if (gVariables)
		for (long i = 0; i<gVariables->lLength; i++)
		{
			long p = gVariables->lData[i];
			_Variable *v = LocateVar (p);
			if (v->IsIndependent())
				l.Insert ((BaseRef)p);
			else
				l2.Insert ((BaseRef)p);
		}
	// additionally, check to see if there is any implicit dependence on the global variables yet unseen
	if (dVariables)
		for (long i = 0; i<dVariables->lLength; i+=2)
		{	
			_SimpleList temp;
			{
				_AVLList  al (&temp);
				_Variable *v = LocateVar (dVariables->lData[i]);			
				v->ScanForVariables(al, true);
				al.ReorderList();
			}
			// see if any of them are global
			for (long j=0; j<temp.lLength;j++)
			{
				long p = temp.lData[j];
				_Variable * v = LocateVar(p);
				if (v->IsGlobal()) // good sign!
				{
					if (v->IsIndependent())
						l.Insert ((BaseRef)p);
					else
						l2.Insert ((BaseRef)p);
					
				}
			}
		}
}

//__________________________________________________________________________________

_Operation::_Operation 	(void)  {numberOfTerms = 0; theData = -1;theNumber = nil;}

//__________________________________________________________________________________

void	_Operation::Initialize(void)
{
	numberOfTerms = 0; theData = -1;theNumber = nil;
}

//__________________________________________________________________________________
BaseRef _Operation::makeDynamic (void)
{
	_Operation * res = new _Operation;
	checkPointer(res);
	//memcpy ((char*)res, (char*)this, sizeof (_Operation));
	res->Duplicate(this);
	return res;
}

//__________________________________________________________________________________
void	_Operation::Duplicate(BaseRef r)
{
	_Operation * o = (_Operation*)r;
	numberOfTerms  = o->numberOfTerms; 
	theData 	   = o->theData;
	theNumber 	   = o->theNumber;
	opCode 		   = o->opCode;
	if (theNumber)
		theNumber->nInstances++;
}


//__________________________________________________________________________________
BaseRef _Operation::toStr(void)
{
	_String res, *dump = nil;
	if (theData!=-1)
	{
		dump = (_String*)((_Variable*)LocateVar(theData))->toStr();
		res = _String("Variable ")& *dump;
	}
	else
	if (theNumber)
	{
		dump = (_String*)theNumber->toStr();
		res = _String("Constant ")& *dump;
	}
	else
		res = _String("Operation ")&*(_String*)BuiltInFunctions(opCode);
	
	if(dump) 
	{
		DeleteObject (dump);
	}
	return res.makeDynamic();
		
}
//__________________________________________________________________________________
_Operation::_Operation 	(_String& opc, long opNo = 2)
// construct the operation by its symbol and, if relevant -
// number of operands
{
	if(opNo>=0)
		opCode = BuiltInFunctions.BinaryFind(&opc);
	else
		opCode = -opNo-1;
		
	if (opCode<0)
	{
		_String errMsg ("Operation: ");
		errMsg = errMsg & opc &" is not defined.";
		WarnError(errMsg);	
	}
	numberOfTerms = opNo;
	theData		  = -1;
	theNumber	  = nil;
}
//__________________________________________________________________________________
_Operation::_Operation 	(_PMathObj theObj)
// construct the operation by its symbol and, if relevant -
// number of operands
{
	numberOfTerms = 0;
	theData		  = -1;
	opCode		  = -1;
	theNumber	  = theObj;
}
//__________________________________________________________________________________
_Operation::_Operation 	(bool isVar, _String& stuff, bool isG, _VariableContainer* theParent)
{
	if (isVar) // creating a variable
	{
		long f;
		_String theS (stuff);
		if (theParent/*&&(!isG)*/) // 20070620: SLKP the commenting may break default behavior!
		{
			f = LocateVarByName(theS);
			
			if (f>=0 && !FetchVar(f)->IsGlobal())
				f = -1;

			if (f<0)
			{
				/*f = stuff.Find('.');
				do 
				{
					if (f!=-1)
					{
						f = stuff.Find(".",f+1,-1); // skip scope levels as needed
						continue;
					}
					else
						break;
				}
				while ((theParent = theParent->theParent));
				if (theParent)*/
					theS = (*theParent->theName)&"."&theS;
			}
			/*StringToConsole(stuff);
			NLToConsole();
			StringToConsole(theS);
			NLToConsole();*/					
		} 
			
		f = LocateVarByName(theS);
		
		if (f<0)
		{
			_Variable v (theS, isG);
			f = v.theIndex;
		}
		else
			f = variableNames.GetXtra(f);			
			
		theData = f;
		theNumber = nil;
		numberOfTerms = 0;
	}
	else
	{
		numberOfTerms = 0;
		theNumber = new _Constant (stuff);
		theData = -1;
	}
	opCode = -1;

}


//__________________________________________________________________________________

_Operation::~_Operation	(void) 
{
	
	if (theNumber)
		DeleteObject (theNumber);
}
//__________________________________________________________________________________
bool _Operation::IsAVariable(bool deep)
{ 
	if (theData==-1)
	{
		if (deep&&theNumber)
			return theNumber->IsVariable();
		return false;
	}
	return true;
		
}

//__________________________________________________________________________________
bool _Operation::IsConstant (void)
{ 
	if (theData==-1)
	{
		if (theNumber)
			return theNumber->IsVariable();
		return true;
	}
	return LocateVar(GetAVariable())->IsConstant();
		
}

//__________________________________________________________________________________

bool		_Operation::EqualOp (_Operation* otherOp)
{
	if (theNumber)
	{
		if (otherOp->theNumber)
		{
			long oc = theNumber->ObjectClass();
			
			if ((oc == NUMBER) && (oc==otherOp->theNumber->ObjectClass()))
				return CheckEqual (theNumber->Value(), otherOp->theNumber->Value());
		}
		return false;
	}
	else
		if (theData==-1)
		{
			if (numberOfTerms<0)
				return numberOfTerms == otherOp->numberOfTerms;
			else
				return opCode == otherOp->opCode;
		}
		else
		{
			return theData == otherOp->theData;
		}
			
	return false;		
}

//__________________________________________________________________________________
	
bool		_Operation::Execute (_Stack& theScrap, _VariableContainer* nameSpace)
{
		if (theNumber)
		{
			theScrap.Push(theNumber);
			return true;
		} 
		if (theData>-1)
		{
			theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[theData])->Compute());
			return true;
		}
		if (theData<-2)
		{
			theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[-theData-3])->GetValue());
			return true;
		}
		if (numberOfTerms<0) // execute a user-defined function
		{
		
			{
				long functionID = -numberOfTerms-1;
				if (functionID >= batchLanguageFunctionParameters.lLength)
				{
					WarnError ("Attempted to call an undefined user function.");
					return false;					
				}
				long numb = batchLanguageFunctionParameters(functionID);
				if (theScrap.StackDepth()<numb) 
				{
					_String errMsg ("User-defined function:");
					errMsg = errMsg&*(_String*)batchLanguageFunctionNames(functionID)
							&" needs "&_String(numb)& " parameters, but only "&_String(theScrap.StackDepth())&" were supplied.";
					WarnError (errMsg);
					return false;
				}
				
				_List 		displacedVars, 
						   *funcVarList = (_List*)batchLanguageFunctionParameterLists(functionID), 
						    displacedValues,
							referenceArgs;
						    
				_SimpleList existingIVars,
							existingDVars, 
							displacedReferences;
				
				_String		*argNameString;
				long		i;
				bool		purgeFlas = false;
				
				for (long k=numb-1; k>=0; k--)
				{
					bool			isRefVar = false;
					argNameString = (_String*)(*funcVarList)(k);
					
					if (argNameString->sData[argNameString->sLength-1]=='&')
					{
						argNameString->Trim(0,argNameString->sLength-2);
						isRefVar  = true;
						purgeFlas = true;
					}
					
					long 	  f = LocateVarByName (*argNameString);
					
					_PMathObj nthterm = theScrap.Pop();
					
					if (isRefVar && nthterm->ObjectClass()!=STRING)
					{
						_FString * type = (_FString*)nthterm->Type();
						_String errMsg = _String ("User-defined function '")
										 &*(_String*)batchLanguageFunctionNames(-numberOfTerms-1)
								         &"' expected a string for the reference variable '" 
										 & *argNameString 
										 &"' but was passed a " & *type->theString 
										 & " with the value of " & _String((_String*)nthterm->toStr());
						
						DeleteObject (type);
						WarnError    (errMsg);
						return false;
					}
					
					if (f<0) // not an existing var
					{
						_Variable newV (*argNameString);
						f =  LocateVarByName (*argNameString);
					}
					_Variable* theV = FetchVar (f);
					
					if (!isRefVar)
					{
						if (theV->IsIndependent())
							// if the variable exists and is independent then 
							// simply swap the value of the var, otherwise
							// duplicate the entire variable
						{
							if (!theV->varValue)
								theV->Compute();
							displacedValues<<theV->varValue;
							theV->varFlags |= HY_VARIABLE_CHANGED;
							/*if (*theV->GetName() == _String("gbdd1"))
							{
								 printf ("Setting argument value of %s to %s\n", 
								 theV->GetName()->sData, 
								 _String((_String*)nthterm->toStr()).sData
							 );
							}*/
							
							theV->varValue = nthterm;
							existingIVars<<theV->GetAVariable();
						}
						else
						{
							_Variable newV (*argNameString);
							newV.SetValue(nthterm);
							existingDVars<<theV->GetAVariable();
							displacedVars<<theV;
						    variablePtrs.Replace (theV->GetAVariable(),(_PMathObj)newV.makeDynamic());
							DeleteObject (nthterm);
						}
					}
					else
					{
						referenceArgs<< variableNames.Retrieve(f);							
						displacedReferences<<theV->GetAVariable();
						 
						_String * refArgName = ((_FString*)nthterm)->theString;
						
						if (nameSpace)
							*refArgName = AppendContainerName (*refArgName, nameSpace);

						i = LocateVarByName (*refArgName);
						
						if (i<0)
						{
							_Variable newAV (*refArgName);
							i =  LocateVarByName (*refArgName);
						}
						variableNames.SetXtra (f, variableNames.GetXtra (i));
						*argNameString = *argNameString&'&';
						DeleteObject (nthterm);
					}
				}
				
				if (purgeFlas)	
					((_ExecutionList*)batchLanguageFunctions(opCode))->ResetFormulae();
					
				_ExecutionList * functionBody = ((_ExecutionList*)batchLanguageFunctions(opCode));
				if (currentExecutionList && currentExecutionList->stdinRedirect)
				{
					functionBody -> stdinRedirect    = currentExecutionList->stdinRedirect;
					functionBody -> stdinRedirectAux = currentExecutionList->stdinRedirectAux;
				}
				_PMathObj ret = functionBody->Execute();
				
				functionBody -> stdinRedirect    = nil;
				functionBody -> stdinRedirectAux = nil;
				
				if (terminateExecution)
				{
					theScrap.Push (new _Constant (0.0));
					return true;
				}
				
				if (ret) 
					theScrap.Push (ret);
				
				for (long di = 0; di < referenceArgs.lLength; di++)
					variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(di)),displacedReferences.lData[di]);				
					
				//for (i= referenceArgs.lLength-1; i>=0; i--)
				//	variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(i)),displacedReferences.lData[i]);				
				
				for (long dv = 0; dv < displacedVars.lLength; dv++)
					variablePtrs.Replace (existingDVars.lData[dv],(_PMathObj)displacedVars(dv));				

				//for (i= displacedVars.lLength-1; i>=0; i--)
				//	variablePtrs.Replace (existingDVars(i),(_PMathObj)displacedVars(i));
					
				for (long dv2 = 0; dv2 < displacedValues.lLength; dv2++)
				{
					_Variable* theV = LocateVar (existingIVars.lData[dv2]);
					DeleteObject(theV->varValue);
					theV->varValue = ((_PMathObj)displacedValues(dv2));
				}
					
				//for (i= displacedValues.lLength-1; i>=0; i--)
				//{
				//	_Variable* theV = LocateVar (existingIVars(i));
				//	DeleteObject(theV->varValue);
				//	theV->varValue = ((_PMathObj)displacedValues(i));
				//}			

		
			}		
			return true;
				
		}

		if (theScrap.theStack.lLength<numberOfTerms) 
		{
			_String *s,errMsg;
			s = (_String*)toStr();
			errMsg = *s;
			DeleteObject(s);
			errMsg = errMsg&
					" needs "&_String(numberOfTerms)& " arguments. Only "&_String(theScrap.StackDepth())&" were given.";
			WarnError (errMsg);
			return nil;
		}
		
		_PMathObj term1, term2 = nil, term3 = nil, temp;

		
		if (numberOfTerms >= 3)
		{
			long sL = theScrap.theStack.lLength-1;
			term3 = (_PMathObj)theScrap.theStack.lData[sL--];
			term2 = (_PMathObj)theScrap.theStack.lData[sL--];		
			term1 = (_PMathObj)theScrap.theStack.lData[sL];
			theScrap.theStack.lLength = sL;		
			temp = term1->Execute (opCode, term2, term3);
			DeleteObject (term1);
			DeleteObject (term2);
			DeleteObject (term3);
		}
		else
			if (numberOfTerms == 2)
			{
				long sL = theScrap.theStack.lLength-1;
				term2 = (_PMathObj)theScrap.theStack.lData[sL--];		
				term1 = (_PMathObj)theScrap.theStack.lData[sL];
				theScrap.theStack.lLength = sL;		
				temp = term1->Execute (opCode, term2, nil);
				DeleteObject (term1);
				DeleteObject (term2);
			}
			else
			{
				term1 = (_PMathObj)theScrap.theStack.lData[--theScrap.theStack.lLength];			
				temp = term1->Execute (opCode, nil, nil);
				DeleteObject (term1);
			}

		
		if (temp)
		{
			theScrap.theStack.Place(temp);
			return true;
		}
		else
			return false;
		
}

//__________________________________________________________________________________
	
void		_Operation::StackDepth (long& depth)
{
	if (theNumber || theData > -1 || theData < -2)
	{
		depth++;
		return;
	}
	
	if (numberOfTerms<0) // execute a user-defined function
	{
		depth -= batchLanguageFunctionParameters(-numberOfTerms-1) - 1;
		return;
	}
	depth -= numberOfTerms - 1;
}

	
//__________________________________________________________________________________

bool		_Operation::ExecutePolynomial (_Stack& theScrap)
{
	if (theData<=-2 || numberOfTerms < 0)
		return false;

	_Polynomial*p = nil;
	if (theNumber)
		p= (_Polynomial*)checkPointer(new _Polynomial(theNumber->Value()));
	
	if (theData>-1)
		p= (_Polynomial*)checkPointer(new _Polynomial(*LocateVar(theData>-1?theData:-theData-2)));
	
	if (p)
	{
		theScrap.Push(p);
		DeleteObject(p);
		return true;
	}
		
	if (theScrap.StackDepth()<numberOfTerms) 
	{
		_String s((_String*)toStr());
		WarnError (s&
				   " needs "&_String(numberOfTerms)& " arguments. Only "&_String(theScrap.StackDepth())&" were given.");
		return 	  false;
	}
	
	_PMathObj term1, 
			  term2 = nil, 
			  temp;

	bool	  opResult = true;
	
	if (numberOfTerms == 2)
		term2 = theScrap.Pop();

	term1 = theScrap.Pop();
	temp  = term1->Execute (opCode, term2, nil);
	DeleteObject (term1);
	
	if (temp)
	{
		theScrap.Push (temp);
		DeleteObject (temp);
	}
	else
		opResult = false;
	
	if (term2)
		DeleteObject (term2);

	return opResult;
}

//__________________________________________________________________________________
	
void	CompileListOfUserExpressions (_SimpleList& varRefs,_List& rec, bool doAll)
{
	rec.Clear();
	if (varRefs.lLength == 0) return;
	
	long i;
	_SimpleList startVars;
	_VariableContainer*  firstVar = (_VariableContainer*)LocateVar(varRefs.lData[0]);
	
	firstVar->ScanAndAttachVariables();
	
	{
		_AVLList sA (&startVars);
		if (doAll)
		{
			
			firstVar->ScanForVariables (sA,sA);
			firstVar->ScanForGVariables (sA,sA);
		}	
		
		firstVar->ScanForDVariables (sA,sA);
		sA.ReorderList ();
	}
	
	if (!doAll)
	{
		for (i=startVars.lLength-1;i>=0;i--)
		{
			if (firstVar->IsModelVar(i))
				startVars.Delete(i);
		}
	}
	
	for (i=0;i<startVars.lLength;i++)
	{
		_String thisName (LocateVar(startVars.lData[i])->GetName()->Cut
						 (LocateVar(startVars.lData[i])->GetName()->FindBackwards('.',0,-1),-1));
		rec && &thisName;
	}
	
	for (i=varRefs.lLength-1;i>=1;i--)
	{
		firstVar = (_VariableContainer*)LocateVar(varRefs.lData[i]);
		firstVar->ScanAndAttachVariables();
		firstVar->MatchParametersToList (rec,doAll);
	}

	for (i=rec.lLength-1;i>=0;i--)
	{
		_String* thisLine = ((_String*)rec(i));
		thisLine->Trim(1,-1);
		if (doAll)
			if (LocateVarByName(*thisLine)<0)
				*thisLine = _String('!')&*thisLine;
	}
	
}
		
//__________________________________________________________
void  FindUnusedObjectName (_String& prefix, _String& partName, _List& names, bool sorted)
{
	if (partName.sLength==0)
		partName = prefix;

	_String tryName (partName);
	long    k = 1;

	if (sorted)
		while (names.BinaryFind(&tryName)>=0)
		{
			k++;
			tryName = partName&k;
		}
	else
		while (names.Find(&tryName)>=0)
		{
			k++;
			tryName = partName&k;
		}
		
	partName = tryName;
}

//__________________________________________________________
void  FindUnusedObjectName (_String& prefix, _String& partName, _AVLListX& names, bool)
{
	if (partName.sLength==0)
		partName = prefix;

	_String tryName (partName);
	long    k = 1;

	while (names.Find(&tryName)>=0)
	{
		k++;
		tryName = partName&k;
	}
	
	partName = tryName;
}

//__________________________________________________________
void  FinishDeferredSF (void)
{
	if (deferSetFormula->lLength)
	{
		SortLists (deferSetFormula, &deferIsConstant);
		_SimpleList tcache;
		long		iv,	
					i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());
					
		for (; i >= 0; i = variableNames.Traverser (tcache,iv))
		{
			_Variable* theV = FetchVar(i);
			if (theV->IsContainer())
				((_VariableContainer*)theV)->SetMDependance (*deferSetFormula);
		}		
		
		for (long j = 0; j<likeFuncList.lLength; j++)
			if (((_String*)likeFuncNamesList(j))->sLength)
			{
				_LikelihoodFunction * lf = (_LikelihoodFunction*)likeFuncList(j);
				for (long k = 0; k < deferSetFormula->lLength; k++)
					lf->UpdateIndependent(deferSetFormula->lData[k],deferIsConstant.lData[k]);
			}
	}
	DeleteObject (deferSetFormula);
	deferSetFormula = nil;
	deferIsConstant.Clear();
}

