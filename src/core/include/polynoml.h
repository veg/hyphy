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

#ifndef		__POLY__
#define		__POLY__

#include "parser.h"

#define	  POLY_DATA_INCREMENT 10
#define	  GLOBAL_VARIABLE   1
#define	  CATEGORY_VARIABLE 2
#define	  RANDOM_VARIABLE   3


//__________________________________________________________________________________

class _PolynomialData : public BaseObj
{

public:

	_PolynomialData (void);
	_PolynomialData (long);
	_PolynomialData (_PolynomialData&);
	_PolynomialData (long,long, _Parameter*);

	virtual	~_PolynomialData ();

	virtual	BaseObj*	makeDynamic(void);
	virtual	void		Duplicate  (BaseRef);

	inline	_Parameter* 		GetCoeff (void) {
		return theCoeff;
	}
	inline	_Parameter& 		GetCoeff (long index) {
		return theCoeff[index];
	}

	long	*   		GetTerm (long);
	long		  		GetNoTerms (void) {
		return actTerms;
	}
	void				AddTerm (long*, _Parameter);
	void				AddTerm (long*, _Parameter, long*, long);
	void				AddTerm (_Parameter);
	void				WriteTerm (long*,long);
	void				DeleteTerm (long);
	bool				IsFirstANumber (void);
	inline long			NumberOfTerms (void) {
		return actTerms;
	}
	long				SumOfPowers (long);
	long				WeightedSumOfPowers (long,_Parameter*);

	// temp!

	bool				checkMe (void);

	friend class _Polynomial;

	void				MultiplyTerms (long*, long*, long*);
	void				RaiseTerm	  (long*, long);
	static	_Parameter	BinaryRaise	  (_Parameter, long);
	static	void		RearrangeTerm (long*, long*, long*,long);
	char				CompareTerms  (long*, long*);
	char				CompareTerms  (long*, long*, long*, long);
	char				CompareTerms  (long*, long*, long*, long*, long, long);
	long				FindTerm	  (long*, long*, long start = 0);
	void				ResortTerms	  (long*);
	void				ChopTerms	  (void);
	bool				checkTerm  	  (_Parameter, long);


protected:

	_Parameter* 	theCoeff;
	long* 			thePowers;
	long   			numberVars, actTerms, allocTerms;

};

//__________________________________________________________________________________

class _Polynomial : public _MathObject
{

public:

	_Polynomial 			(void);
	_Polynomial 			(_SimpleList&);
	_Polynomial 			(_Polynomial&);
	_Polynomial 			(_Parameter);
	_Polynomial 			(_Variable&);
	virtual					~_Polynomial ();
	virtual 				_MathObject* Execute (long opCode, _MathObject* p = nil , _MathObject* p2 = nil);   // execute this operation with the list of Args

	virtual	BaseObj*		makeDynamic(void);
	virtual	void			Duplicate  (BaseRef);

	virtual	_MathObject* 	Add 				(_MathObject*);
	virtual	_MathObject* 	Plus 				(_MathObject*, bool subtract = false);
	virtual	_MathObject* 	Sub 				(_MathObject*);
	virtual	_MathObject* 	Raise 				(_MathObject*);
	virtual	_MathObject* 	Minus 				(void);
	virtual	_MathObject* 	Mult 				(_MathObject*);
	virtual	_MathObject* 	Compute 			(void);
	virtual	bool			Equal				(_MathObject*);
	_Parameter				ComputePolynomial	(void);

	_Parameter 				ComputeP 			(_Parameter* , _Parameter* , long , long, long*, long*);
	_MathObject*			IsANumber			(bool = false);
	virtual	 bool 			IsObjectEmpty 		(void);

	virtual	long			ObjectClass (void) {
		return POLYNOMIAL;
	}
	virtual	_Parameter		Value (void) {
		return ComputePolynomial();
	}

	virtual BaseObj* 		toStr (void);
	void			 		CheckTerm(void);

	virtual void	 		toFileStr (FILE*);

	long					GetNoVariables(void) {
		return variableIndex.countitems();
	}
	_PolynomialData*		GetTheTerms(void) {
		return theTerms;
	}
	void					SetTheTerms(_PolynomialData* td) {
		theTerms = td;
	}
	void					SetCLists(_SimpleList& c1,_SimpleList& c2) {
		compList1.Duplicate(&c1);
		compList2.Duplicate(&c2);
	}
	virtual	void			ScanForVariables
	(_AVLList &l, bool globals = false);
	virtual	bool			HasChanged (void);
	friend	void			ResetPolynomialCheck
	(_Polynomial*);
	long					ComputationalSize (void) {
		return compList1.countitems();
	}
	bool					IsMaxElement 	(_Parameter);
	void		 			Convert2ComputationForm
	(_SimpleList *c1 = nil, _SimpleList *c2 = nil, _SimpleList* termsToInclude = nil);
	void					RankTerms 		(_SimpleList*);
protected:

	void		 		DropSmallTerms(void);
	void		 		Convert2OperationForm(void);

	_SimpleList	 		variableIndex,
						compList1,
						compList2;

	_PolynomialData	*   theTerms;


};

extern _Parameter dropPrecision, topPolyCap, dropTerms, enforcePolyCap,
	   maximumPolyTermsPerVariable, maxPolynomialExpIterates,polynomialExpPrecision;
void	SetPolyTermCap (long);
#endif