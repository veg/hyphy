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

#ifndef		__PARSER__
#define		__PARSER__

#include "baseobj.h"
#include "hy_lists.h"
#include "hy_strings.h"
#include "errorfns.h"
#include "stdio.h"
#include "classes.h"


#define  HY_UNDEFINED 		0x00
#define  NUMBER 			0x01
#define	 MATRIX 			0x04
#define	 CONTAINER			0x08
#define	 TREE_NODE			0x10
#define	 TREE				0x20
#define	 STRING				0x40
#define	 ASSOCIATIVE_LIST 	0x80
#define	 TOPOLOGY			0x100
#define	 POLYNOMIAL			0x200
#define  HY_ANY_OBJECT		0xFFFF

#define	 DEFAULTLOWERBOUND 	-1e26 
#define	 DEFAULTUPPERBOUND	1e26

#define	 GLOBAL_VARIABLE   1
#define	 CATEGORY_VARIABLE 2
#define	 RANDOM_VARIABLE   3

#define	 HY_VARIABLE_GLOBAL 	0x01
#define	 HY_VARIABLE_CHANGED	0x02
#define	 HY_DEP_V_COMPUTED	    0x04
#define	 HY_DEP_V_INSPECTED	    0x08
#define	 HY_DEP_V_INSPECTED_CLR	0xF7
#define	 HY_DEP_V_MODIFIED	    0x10
#define	 HY_DEP_V_MODIFIED_CATS	0x20
#define	 HY_VC_NO_CHECK			0x40 // do not check this variable container in 
									 // NeedToExponentiate
#define	 HY_VARIABLE_NOTSET 	0x80 
#define	 HY_VARIABLE_SET		0x7F

#define	 HY_VC_CLR_NO_CHECK		0xBF 

#define  HY_DEP_CLEAR_MASK		0xC7

#define  HY_NO_MODEL			(-1)


// START OPCODES

#define	 HY_OP_CODE_NOT				0								// !
#define	 HY_OP_CODE_NEQ				(1+HY_OP_CODE_NOT)				// !=
#define	 HY_OP_CODE_IDIV			(1+HY_OP_CODE_NEQ)				// $
#define	 HY_OP_CODE_MOD				(1+HY_OP_CODE_IDIV)				// %
#define	 HY_OP_CODE_AND				(1+HY_OP_CODE_MOD)				// &&
#define	 HY_OP_CODE_MUL				(1+HY_OP_CODE_AND)				// *
#define	 HY_OP_CODE_ADD				(1+HY_OP_CODE_MUL)				// +
#define	 HY_OP_CODE_SUB				(1+HY_OP_CODE_ADD)			// -
#define	 HY_OP_CODE_DIV				(1+HY_OP_CODE_SUB)				// /
#define	 HY_OP_CODE_LESS			(1+HY_OP_CODE_DIV)  // <

#define	 HY_OP_CODE_LEQ				(1+HY_OP_CODE_LESS) // <=
#define	 HY_OP_CODE_EQ				(1+HY_OP_CODE_LEQ) // ==
#define	 HY_OP_CODE_GREATER			(1+HY_OP_CODE_EQ) // >
#define	 HY_OP_CODE_GEQ				(1+HY_OP_CODE_GREATER) // >=
#define	 HY_OP_CODE_ABS				(1+HY_OP_CODE_GEQ) // Abs
#define	 HY_OP_CODE_ARCTAN			(1+HY_OP_CODE_ABS) // Arctan
#define	 HY_OP_CODE_BETA			(1+HY_OP_CODE_ARCTAN) // Beta
#define	 HY_OP_CODE_BRANCHCOUNT		(1+HY_OP_CODE_BETA) // BranchCount
#define	 HY_OP_CODE_BRANCHLENGTH	(1+HY_OP_CODE_BRANCHCOUNT) // BranchLength
#define	 HY_OP_CODE_BRANCHNAME		(1+HY_OP_CODE_BRANCHLENGTH) // BranchName

#define	 HY_OP_CODE_CCHI2			(1+HY_OP_CODE_BRANCHNAME) // CChi2
#define	 HY_OP_CODE_CGAMMADIST		(1+HY_OP_CODE_CCHI2) // CGammaDist
#define	 HY_OP_CODE_COLUMNS			(1+HY_OP_CODE_CGAMMADIST) // Columns
#define	 HY_OP_CODE_COS				(1+HY_OP_CODE_COLUMNS) // Cos
#define	 HY_OP_CODE_EIGENSYSTEM		(1+HY_OP_CODE_COS) // Eigensystem
#define	 HY_OP_CODE_ERF				(1+HY_OP_CODE_EIGENSYSTEM) // Erf
#define	 HY_OP_CODE_EVAL			(1+HY_OP_CODE_ERF) // Eval
#define	 HY_OP_CODE_EXP				(1+HY_OP_CODE_EVAL) // Exp
#define	 HY_OP_CODE_FORMAT			(1+HY_OP_CODE_EXP) // Format
#define	 HY_OP_CODE_GAMMA			(1+HY_OP_CODE_FORMAT) // Gamma

#define	 HY_OP_CODE_GAMMADIST		(1+HY_OP_CODE_GAMMA) // GammaDist
#define	 HY_OP_CODE_IBETA			(1+HY_OP_CODE_GAMMADIST) // IBeta
#define	 HY_OP_CODE_IGAMMA			(1+HY_OP_CODE_IBETA) // IGamma
#define	 HY_OP_CODE_INVCHI2			(1+HY_OP_CODE_IGAMMA) // InvChi2
#define	 HY_OP_CODE_INVERSE			(1+HY_OP_CODE_INVCHI2) // Inverse
#define	 HY_OP_CODE_LUDECOMPOSE		(1+HY_OP_CODE_INVERSE) // LUDecompose
#define	 HY_OP_CODE_LUSOLVE			(1+HY_OP_CODE_LUDECOMPOSE) // LUSolve
#define	 HY_OP_CODE_LOG				(1+HY_OP_CODE_LUSOLVE) // Log
#define	 HY_OP_CODE_MACCESS			(1+HY_OP_CODE_LOG) // MAccess
#define	 HY_OP_CODE_MCOORD			(1+HY_OP_CODE_MACCESS) // MCoord

#define	 HY_OP_CODE_MAX				(1+HY_OP_CODE_MCOORD) // Max
#define	 HY_OP_CODE_MIN				(1+HY_OP_CODE_MAX) // Min
#define	 HY_OP_CODE_PSTREESTRING	(1+HY_OP_CODE_MIN) // PSTreeString
#define	 HY_OP_CODE_RANDOM			(1+HY_OP_CODE_PSTREESTRING) // Random
#define	 HY_OP_CODE_REROOTTREE		(1+HY_OP_CODE_RANDOM) // RerootTree
#define	 HY_OP_CODE_ROWS			(1+HY_OP_CODE_REROOTTREE) // Rows
#define	 HY_OP_CODE_SIMPLEX			(1+HY_OP_CODE_ROWS) // Simplex
#define	 HY_OP_CODE_SIN				(1+HY_OP_CODE_SIMPLEX) // Sin
#define	 HY_OP_CODE_SQRT			(1+HY_OP_CODE_SIN) // Sqrt
#define	 HY_OP_CODE_TEXTREESTRING	(1+HY_OP_CODE_SQRT) // TEXTreeString

#define	 HY_OP_CODE_TAN				(1+HY_OP_CODE_TEXTREESTRING) // Tan
#define	 HY_OP_CODE_TIME			(1+HY_OP_CODE_TAN) // Time
#define	 HY_OP_CODE_TIPCOUNT		(1+HY_OP_CODE_TIME) // TipCount
#define	 HY_OP_CODE_TIPNAME			(1+HY_OP_CODE_TIPCOUNT) // TipName
#define	 HY_OP_CODE_TRANSPOSE		(1+HY_OP_CODE_TIPNAME) // Transpose
#define	 HY_OP_CODE_TYPE			(1+HY_OP_CODE_TRANSPOSE) // Type
#define	 HY_OP_CODE_ZCDF			(1+HY_OP_CODE_TYPE) // ZCDF
#define	 HY_OP_CODE_POWER			(1+HY_OP_CODE_ZCDF) // ^
#define	 HY_OP_CODE_OR				(1+HY_OP_CODE_POWER) // ||

// END OPCODES

class   _Variable;
class   _VariableContainer;

class	_MathObject : public BaseObj{ //abstract math operations class

	public:

	virtual	_MathObject* Add 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Sub 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Minus 		(void) 		   	   {return nil;}
	virtual	_MathObject* Mult 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Div 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* lDiv 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* longDiv 	(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Raise 		(_MathObject*) 	   {return nil;}
	virtual void		 Assign 	(_MathObject*) 	   {}
	virtual	bool		 Equal 		(_MathObject*) 	   {return false;}
	virtual	_MathObject* Abs 		(void) 			   {return nil;}
	virtual	_MathObject* Sin 		(void) 			   {return nil;}
	virtual	_MathObject* Cos 		(void) 			   {return nil;}
	virtual	_MathObject* Tan 		(void) 			   {return nil;}
	virtual	_MathObject* Exp 		(void) 			   {return nil;}
	virtual	_MathObject* Log 		(void) 			   {return nil;}
	virtual	_MathObject* Sqrt 		(void) 			   {return nil;}
	virtual	_MathObject* Gamma 		(void) 			   {return nil;}
	virtual	_MathObject* Erf 		(void) 			   {return nil;}
	virtual _MathObject* LnGamma	(void)			   {return nil;}	// <-- added by afyp, February 7, 2007
	virtual	_MathObject* Beta  		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* IGamma		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* CChi2		(_MathObject*) 	   {return nil;}
	virtual _MathObject* IBeta		(_MathObject*,_MathObject*)
													   {return nil;}
	virtual _MathObject* Simplex	(void)			   {return nil;}
	virtual _MathObject* Min		(_MathObject*)	   {return nil;}
	virtual _MathObject* Max		(_MathObject*)	   {return nil;}
	virtual	_MathObject* InvChi2	(_MathObject*) 	   {return nil;}
	virtual	_MathObject* ZCDF		(void)		 	   {return nil;}
	virtual	_MathObject* Time 		(void) 		   	   {return nil;}
	virtual	_MathObject* Arctan 	(void) 		   	   {return nil;}
	virtual	_MathObject* Less 		(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Random 	(_MathObject*) 	   {return nil;}
	virtual	_MathObject* Greater	(_MathObject*)     {return nil;}
	virtual	_MathObject* LessEq 	(_MathObject*) 	   {return nil;}
	virtual	_MathObject* GreaterEq 	(_MathObject*)     {return nil;}
	virtual	_MathObject* AreEqual 	(_MathObject*)     {return nil;}
	virtual	_MathObject* NotEqual 	(_MathObject*)     {return nil;}
	virtual	_MathObject* LAnd	 	(_MathObject*)     {return nil;}
	virtual	_MathObject* LOr	 	(_MathObject*)     {return nil;}
	virtual	_MathObject* GammaDist  (_MathObject*,_MathObject*)  
													   {return nil;}
	virtual	_MathObject* CGammaDist (_MathObject*,_MathObject*) 
													   {return nil;}
	virtual	_MathObject* LNot  		(void)	 		   {return nil;}
	virtual	_MathObject* TipCount  	(void)		   	   {return nil;}
	virtual	_MathObject* BranchCount (void)        	   {return nil;}
	virtual	_MathObject* TipName	 (_MathObject*)	   {return nil;}
	virtual	_MathObject* BranchName	 (_MathObject*)	   {return nil;}
	virtual	_MathObject* BranchLength(_MathObject*)	   {return nil;}
	virtual	_MathObject* RerootTree	 (_MathObject*)	   {return nil;}
	virtual	_MathObject* TEXTreeString(_MathObject*)
												   	   {return nil;}
	virtual	_MathObject* Type						   (void);
	virtual	_MathObject* PlainTreeString(_MathObject*)
						 						   	   {return nil;}
	virtual	_MathObject* FormatNumberString (_MathObject*,_MathObject*)
												   {return nil;}
	virtual	_Parameter	 Value (void) 			   {return 0.0;}
	virtual	_MathObject* Compute (void) 		   {return nil;}
	virtual	void	     ScanForVariables (_AVLList&,bool = false)			
													{}
	
	virtual	bool		 IsVariable (void)		   {return false;}
	virtual	bool		 IsObjectEmpty (void)	   {return true;}
	virtual	bool		 IsPrintable (void)		   {return false;}
	
	virtual	bool		 IsDefined	(_String&);  // is this operation defined for the type
											
	virtual	bool		 IsIndependent (void) 		{ return true; } 
	virtual	long		 ObjectClass (void) 		{ return HY_UNDEFINED; } 
			// returns a unique ID for this object
			// 0 - undefined
			// 1 - number
			// 4 - matrix
											
	virtual _MathObject* Execute (long opCode, _MathObject* p = nil , _MathObject* p2 = nil);   // execute this operation with the list of Args
	virtual	bool		 HasChanged (void) { return false; }

	virtual   bool 	  IsConstant	(void)
						{
							return true;
						}
};

typedef	_MathObject* _PMathObj ;

// pointer to a math object

//__________________________________________________________________________________

class	_Stack { //computational stack

	friend class _Formula;
	friend class _Operation;
	friend class _Variable;

	public:
	
	_Stack (void);
	~_Stack (void);
	
	bool	  Push (_PMathObj); 	// push object onto the stack
	_PMathObj Pop (bool del = true); 			// pop object from the top of the stack
	long	  StackDepth (void); 	// returns the depth of the stack
	void	  Reset (void); 		// clear the stack
	
	virtual	  void	  Initialize (void); 
	virtual	  void	  Duplicate (BaseRef); 
	
	protected:
	
		_List  theStack;
};
	
//__________________________________________________________________________________

class _Constant : public _MathObject { // a numerical constant
	
	public:

	_Constant (_Parameter);
	_Constant (_String&);
	_Constant (void);
	~_Constant (void) {}
	
	virtual	_PMathObj Add 			(_PMathObj);
	virtual	_PMathObj Sub 			(_PMathObj);
	virtual	_PMathObj Minus	 		(void) ;
	virtual	_PMathObj Mult 			(_PMathObj);
	virtual	_PMathObj Div 			(_PMathObj);
	virtual	_PMathObj lDiv 			(_PMathObj);
	virtual	_PMathObj longDiv 		(_PMathObj);
	virtual	_PMathObj Raise 		(_PMathObj);
	virtual void	  Assign 		(_PMathObj);
	virtual	bool	  Equal 		(_PMathObj);
	virtual	_PMathObj Abs 			(void);
	virtual	_PMathObj Sin 			(void);
	virtual	_PMathObj Cos 			(void);
	virtual	_PMathObj Tan 			(void);
	virtual	_PMathObj Exp 			(void);
	virtual	_PMathObj Log 			(void);
	virtual	_PMathObj Sqrt 			(void);
	virtual	_PMathObj Time 			(void);
	virtual	_PMathObj Arctan 		(void);
	virtual	_PMathObj Gamma  		(void);
	virtual _PMathObj LnGamma		(void);			/* <- added by afyp, February 8, 2007 */
	virtual	_PMathObj Beta   		(_PMathObj);
	virtual	_PMathObj Min    		(_PMathObj);
	virtual	_PMathObj Max    		(_PMathObj);
	virtual	_PMathObj GammaDist   	(_PMathObj,_PMathObj);
	virtual	_PMathObj CGammaDist    (_PMathObj,_PMathObj);
	virtual	_PMathObj IBeta	        (_PMathObj,_PMathObj);
	virtual	_PMathObj IGamma   		(_PMathObj);
	virtual	_PMathObj CChi2   		(_PMathObj);
	virtual	_PMathObj InvChi2   	(_PMathObj);
	virtual	_PMathObj Erf   		(void);
	virtual	_PMathObj ZCDF   		(void);
	virtual	_PMathObj Less    		(_PMathObj);
	virtual	_PMathObj Greater 		(_PMathObj);
	virtual	_PMathObj LessEq 		(_PMathObj);
	virtual	_PMathObj GreaterEq 	(_PMathObj);
	virtual	_PMathObj AreEqual 		(_PMathObj);
	virtual	_PMathObj NotEqual 		(_PMathObj);
	virtual	_PMathObj LAnd	 		(_PMathObj);
	virtual	_PMathObj LOr	 		(_PMathObj);
	virtual	_PMathObj LNot	 		();
	virtual	_PMathObj Random 		(_PMathObj);
	virtual	_Parameter
					  Value 		(void);
	virtual	_PMathObj FormatNumberString 
									(_PMathObj,_PMathObj);
	virtual	_PMathObj Compute 		(void) 
									{return this;};
	
	virtual	  void	  Initialize 			(void); 
	virtual	  void	  Duplicate 			(BaseRef); 
	virtual	  BaseRef makeDynamic			(void);
	virtual	  BaseRef toStr 				(void);
	virtual	  long 	  ObjectClass 			(void) 
						{ return NUMBER;}
	virtual	  void	  SetValue 				(_Parameter pl) 
						{theValue = pl;} 
	
	public:
		
		_Parameter theValue;
		
};

extern	_List BuiltInFunctions;

//__________________________________________________________________________________

class	_Operation : public BaseObj {

	friend class _Formula; 
	friend class _Variable; 
	friend class _VariableContainer; 

	public:

	
	_Operation 	(void);
	_Operation 	(_String&, long); 
							// construct the operation by its symbol and, if relevant -
						    // number of operands
	
	_Operation	(bool, _String&, bool isG = false, _VariableContainer*  = nil); 
							// store a variable or a constant
	_Operation	(_PMathObj); 		
							// store a non-numeric constant
	
	virtual	~_Operation	(void);
	
virtual	  BaseObj*		makeDynamic 		(void);
	
		  bool			Execute 			(_Stack&, _VariableContainer* = nil); //execute this operation
		// see the commend for _Formula::ExecuteFormula for the second argument
virtual   void			StackDepth			(long&);
		  								
		  bool			ExecutePolynomial 	(_Stack&);	
virtual	  BaseObj*		toStr				(void);    //convert the op to string
										
virtual	  void			Initialize			(void); 
virtual	  void			Duplicate 			(BaseRef); 
		  _String&  	GetCode				(void) 
		  					{return (opCode>-1)&&(numberOfTerms>=0)?*(_String*)BuiltInFunctions(opCode):empty;}
		 long&  		TheCode 			(void) 
		 					{return opCode;}
virtual	 bool 			IsAVariable 		(bool = true) ;	// is this object a variable or not?	
virtual	 bool			IsConstant			(void); 		// does this object depend on any independent variables or not?				 	

virtual	 long			UserFunctionID		(void) { return numberOfTerms < 0 ? -numberOfTerms-1 : -1;};	
						// return a non-neg number (function index) if this is a user function,
						// otherwise, return -1

virtual	 long    		GetAVariable 		(void) 		// return the index of the variable
							{return theData>=-2?theData:-theData-3; } 
									
virtual	 void    		SetAVariable 		(long d) 	// return the index of the variable
							{theData=d;} 	
							
virtual  bool			AssignmentVariable	(void)
							{return theData<-2;}
							
virtual	 void   	 	SetTerms 			(long d) 
							{numberOfTerms=d;} 
								
virtual	 _PMathObj 		GetANumber 			(void) 
							{return theNumber;} 
							
virtual	 void   	 	SetNumber 			(_PMathObj d) 
							{theNumber=d;} 
																
		long			GetNoTerms			(void) 
							{return numberOfTerms;}
		long			PrecedenceLevel		(void);

		
virtual bool			EqualOp				(_Operation*);
	
	protected:
	
		long	 	opCode; 		// internal operation code
		long		numberOfTerms,  // 1 - unary, 2 - binary, etc
				 	theData;
		_PMathObj   theNumber;
};
	
	
	
//__________________________________________________________________________________

union		_SimpleFormulaDatum
{
	_Parameter value;
	Ptr		   reference;
};


//__________________________________________________________________________________

class 	_Formula { // a computational formula

	friend class _Variable;
	friend class _VariableContainer; 

	public:
	
			_Formula (void);
			_Formula (_String&,_VariableContainer* theParent=nil,bool errors=true);
			_Formula (_PMathObj, bool isAVar = false);
	
			~_Formula (void);
	
		_PMathObj 	Compute 			(long = 0, _VariableContainer* = nil); 
					// compute the value of the formula
					// 1st argument : execute from this instruction onwards 
					// see the commend for ExecuteFormula for the second argument
	
		bool	 	IsEmpty				(void); // is there anything in the formula
		long	 	NumberOperations 	(void); // how many ops in the formula?
	
friend	long	  	Parse 				(_Formula*, _String&, _VariableContainer* = nil, _Formula* = nil, bool flagErrors = true); // the parser
friend	long	  	ExecuteFormula 		(_Formula*, _Formula*, long, _VariableContainer* = nil); 
										// the execution block for "compiled formulae
	/*  
	 SLKP 20100119: added an execution name space to allow correct scoping of "pass-by-reference" 
					arguments when calling ExecuteAFile within a namespace.
	 
					e.g. in  
	 
					function foo (var&)
				    {
						...
					}
	 
					foo ("varID");
					
					varID may need to be prefixed by a namespace ID.
	 */
	
		bool		CheckFormula 		(void); // check to see if this formula is valid and compute the obj class
	
		_MathObject*ConstructPolynomial (void);
	
virtual	void	  	Initialize 			(void); 
virtual	void	  	Duplicate 			(BaseRef); 
		void	  	DuplicateReference 	(_Formula*); 
virtual	BaseRef		makeDynamic			(void);
virtual BaseRef 	toStr 				(_List*	matchNames = nil, bool = false);
	
virtual	long	  	ObjectClass 		(void);
					
	
virtual	void	  	ScanFForVariables   (_AVLList&l, bool includeGlobals = false, bool includeAll = false, bool includeCateg = true, bool = false);
virtual	void	  	ScanFForType		(_SimpleList&,  int);
		/* SLKP 20100716: 
				A simple utility function to retrieve all variables of a given type
		 */
	
virtual	bool	  	CheckFForDependence (long, bool checkAll = false);
		_List&	  	GetList 			(void) 
						{return theFormula;}	
		
		bool	  	HasChanged 			(bool = false); // does  the formula need recomputing
		bool		EqualFormula		(_Formula*);
		bool	  	IsAConstant 		(void); //  does this formula include variables, or is it just a constant? 
		bool	  	IsConstant 			(void); //  does this formula depend on something other that constants and fixed parameters? 
		bool	  	DependsOnVariable 	(long); 
					/*  
						SLKP 20090315: added a missing utility function
						given a variable index as an argument, returns true if
						the formula depends on a it; false otherwise
					*/
		_Operation*	GetIthTerm			(long); 
					/*  
						SLKP 20090315: added a missing utility function
						given an index (i) as the argument, the function retrieves
						the i-th term of the formula
					*/
		void	  	Clear 				(void);
		_PMathObj 	GetTheMatrix		(void);
	
		bool	  	AmISimple	  		(long& stackDepth, _SimpleList& variableIndex);
		void 	   	ConvertToSimple 	(_SimpleList& variableIndex);
		void 	   	ConvertFromSimple 	(_SimpleList& variableIndex);
		void	   	SimplifyConstants	(void);
	
		_Parameter	ComputeSimple 		(_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues) ;
	
		_Parameter	Newton 				(_Formula&, _Variable*,  _Parameter, _Parameter, _Parameter);
		_Parameter 	Newton 				(_Formula&, _Parameter, _Parameter, _Parameter, _Variable*);
		_Parameter	Newton 				(_Variable*,  _Parameter, _Parameter, _Parameter, _Parameter);
		_Parameter 	Newton 				(_Variable*,_Parameter, _Parameter, _Parameter);

		_Parameter 	Brent 				(_Variable*, _Parameter, _Parameter, _Parameter = 1.e-7, _List* = nil);

		_Parameter	Integral 			(_Variable*,_Parameter, _Parameter, bool inifinite = false);
		_Parameter	MeanIntegral 		(_Variable*,_Parameter, _Parameter, bool inifinite = false);
		_Formula*   Differentiate		(_String, bool = true);
		node<long>*	InternalDifferentiate
										(node<long>*, long,_SimpleList&, _SimpleList&, _Formula&);
										
		bool		InternalSimplify	(node<long>*);
	
		void	    LocalizeFormula 	(_Formula&, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv);
		
	protected:
	
		void		internalToStr 		(_String& result,node<long>*, char opLevel, _List* matchNames, _Operation* = nil);
		void    	ConvertToTree 		(void);
		void		ConvertFromTree 	(void);
		bool		CheckSimpleTerm		(_PMathObj);
		node<long>* DuplicateFormula 	(node<long>*,_Formula&);  
		
		_List  		theFormula;
		_Stack 		theStack;
		node<long>* theTree; // this formula converted to a tree for operation purposes
							 // such as simplification, differentiation and printing.
							 // trees store numbers referencing operations inside
							 // "theFormula"

};
	
//__________________________________________________________________________________
	
class _Variable : public _Constant {

	friend class _Operation;
	
	public:

	_Variable (void);
	_Variable (_String&, bool isG = false); // name
	_Variable (_String&, _String&, bool isG = false); // name and formula
	
	virtual ~_Variable (void);
	
	virtual	  void	  		Initialize (void); 
	virtual	  void	  		Duplicate (BaseRef); 
	virtual	  BaseRef 		makeDynamic(void);
	virtual	  BaseRef 		toStr (void);
	virtual    void	 		toFileStr (FILE*);

	virtual	  void	  		MarkDone (void); 

	virtual		_PMathObj   Compute (void); 	  // compute or return the value
	virtual		bool	    IsVariable (void); //  
	virtual		bool	    IsIndependent (void) 
									{ return (varFormula&&varFormula->theFormula.lLength)?
														false:
														(varValue?varValue->IsIndependent():true);
									}   
	virtual		bool	    IsConstant (void);
				void	    SetValue (_PMathObj, bool = true); // set the value of the variable
				void	    SetNumericValue (_Parameter);
				void	    CheckAndSet (_Parameter, bool = false); 
									// set the value of the variable
									// bool flag is used to indicate that out of bounds values should be rejected
	
				_PMathObj   GetValue (void) {return varValue;} // get the value of the variable
				void	    SetFormula (_Formula&); // set the variable to a new formula
	
	virtual		bool	    HasChanged     (bool = false);
	virtual     void	    PreMarkChanged  ();
	virtual		void	    PostMarkChanged ();
	virtual		bool	    IsGlobal (void) 
											{ return varFlags & HY_VARIABLE_GLOBAL;}
	virtual		bool	    IsCategory (void) 
						 					{ return false;}
	virtual		long	    GetAVariable (void) 
											{ return theIndex;}
	virtual		long	    ObjectClass (void) 
											{ return varValue?varValue->ObjectClass():((varFormula&&varFormula->theFormula.lLength)?varFormula->ObjectClass():1);}
				void	    SetIndex (long i) 
											{theIndex = i;}
				long	    GetIndex (void) 
											{ return theIndex;}
	virtual		void	    ScanForVariables (_AVLList& l, bool globals = false)
						  					{ 
						  						if (varValue) 
						  							varValue->ScanForVariables (l, globals);
						   						if (varFormula && varFormula->theFormula.lLength) 
						   							varFormula->ScanFForVariables(l,globals);
						  					}
						   
	virtual		bool	    IsContainer (void) 
											{ return false;}
	
				void	    SetBounds (_Parameter lb, _Parameter ub) 
											{lowerBound = lb; upperBound = ub;}
				
				_Parameter 	GetLowerBound (void) 
											{ return lowerBound; }
				_Parameter 	GetUpperBound (void) 
											{ return upperBound; }
				
	virtual		void		ClearConstraints 	(void);
	virtual		bool	  	CheckFForDependence (long, bool = false);
				
				_String*	GetName					(void)	
							{ return theName;}
				_String*	GetFormulaString		(void)	
							{return varFormula?(_String*)varFormula->toStr():(_String*)empty.makeDynamic();}
							
	virtual		void		CompileListOfDependents (_SimpleList&);
				
	
	friend		void	    ResetVariables			(void);
	friend		_Variable*  LocateVar 				(long);
	friend		void	    InsertVar 				(_Variable*);
				
	public:
	
		_String*   theName;
		
		_PMathObj  varValue;
		
		long	   theIndex; // index of this variable in the global variable pool

		// the class of this variable - i.e global, local, category or random
		char	   varFlags;
		
		_Parameter lowerBound, 
				   upperBound;
				   // dynamic lower and upper bounds here
				   
		_Formula*  varFormula;
		
};

#include "matrix.h"

//__________________________________________________________________________________

#define	 USE_POINTER_VC

//__________________________________________________________________________________

// this class defines a computational (or storage) class which, as a variable, may contain
// other variables locally. 

class	_VariableContainer: public _Variable {

	friend class _Operation;
	friend class _Variable;
	
	public:
	
	_VariableContainer (void);
	_VariableContainer (_String theName, _String theTmplt = "", _VariableContainer* theP = nil);
	// name, matrix constructor, the parent (if there is one)
	virtual ~_VariableContainer(void);

	void					InitializeVarCont (_String&, _String&, _VariableContainer*, _AVLListXL* = nil);
	
	virtual	    void		MarkDone (void); 

	// variable access/operation functions
	
	virtual		bool	    IsContainer 				(void) 
														{return true;}
					
	virtual		bool	    HasChanged 					(void);
	virtual		bool	    NeedToExponentiate 			(bool = false);
	
				void 	    ScanAndAttachVariables 		(void); 
	
	virtual		void	    ScanForVariables 			(_AVLList&,_AVLList&); 
	virtual		void	    ScanForDVariables 			(_AVLList&,_AVLList&); 
	virtual		void	    ScanForGVariables 			(_AVLList&,_AVLList&); 
	
	virtual		bool	    IsModelVar					(long);
	virtual		bool	    IsConstant					(void);
	virtual		BaseRef	    makeDynamic 				(void);
	virtual		void	    Duplicate   				(BaseRef);
	
	virtual	    BaseRef	    toStr						(void);
	
				bool	    HasLocals 					(void);
	
	virtual		bool		RemoveDependance 			(long);
	virtual		long		SetDependance  				(long);
				bool		SetMDependance 				(_SimpleList&);
	
				void		Clear 						(void);
	virtual		void		ClearConstraints 			(void);
	
				long		CountIndependents 			(void);	
				long		CountAll 					(void);
			
	virtual		_Variable*  GetIthIndependent 			(long);	
	virtual		_Variable*  GetIthDependent 			(long);	
	virtual		_Variable*  GetIthParameter 			(long);	
	
				long		CheckAndAddUserExpression 	(_String&, long startWith = 0);
				void		KillUserExpression 			(long);
	virtual		void		CompileListOfDependents 	(_SimpleList&);
	
				void		MatchParametersToList		(_List&, bool doAll = false, bool indOnly = false);
				_Matrix*	GetModelMatrix 				(void);	
				_Matrix*	GetFreqMatrix 				(void);	
				bool		HasExplicitFormModel		(void);

				long		GetModelIndex 				(void) 
														{ return theModel; }
				long		GetModelDimension			(void);
					/* 20100316 SLKP
						return the dimension of the model; needed to handle the case
						of explicit model exponentials
					 */
	
				void		CopyMatrixParameters 		(_VariableContainer*);
				void		GetListOfModelParameters 	(_List&);
				_String*	GetSaveableListOfUserParameters 		
														(void);
				void		TrimMemory 					(void);
		_VariableContainer* GetTheParent				(void)
														{
															return theParent;
														}

	protected: // data members
			
		_SimpleList		    *iVariables, 
							*dVariables,
							*gVariables;
		
		void			   SortVars (void);

		long				theModel; 	// model template for the container
		_VariableContainer  *theParent; // a higher level container, if there is one.

};					 
		
//__________________________________________________________________________________

class _FString : public _MathObject { // strings encountered in formulas
	
	public:

	_FString (_String&, bool = true);
	_FString (long);
	_FString (_String*);
	_FString (void);
	virtual	 ~_FString ();
//	~_Constant (void);
	
	virtual	BaseRef	  makeDynamic 		(void);
	virtual	void      Duplicate 		(BaseRef);
	virtual	_PMathObj Add 				(_PMathObj);
	virtual	long	  AddOn	 			(_PMathObj);
	virtual	_PMathObj AreEqual 			(_PMathObj);
	virtual	_PMathObj AreEqualCIS 		(_PMathObj);
	virtual	_PMathObj Less	 			(_PMathObj);
	virtual	_PMathObj LessEq	 		(_PMathObj);
	virtual	_PMathObj Greater	 		(_PMathObj);
	virtual	_PMathObj GreaterEq	 		(_PMathObj);
	virtual	_PMathObj NotEqual 			(_PMathObj);
	virtual	_PMathObj RerootTree 		(void);
	virtual	_PMathObj EqualAmb	 		(_PMathObj);
	virtual	_PMathObj EqualRegExp	 	(_PMathObj,bool = false);
	virtual	_PMathObj ReplaceReqExp	 	(_PMathObj);
	virtual	_PMathObj CountGlobalObjects(void);
	virtual _PMathObj FileExists		(void);
	virtual _PMathObj Evaluate			(void);
	virtual	long 	  ObjectClass 		(void) 
					  					{ return STRING;}
	virtual	_PMathObj Compute 			(void) 
					  					{ return this; }
	
	virtual _PMathObj MapStringToVector (_PMathObj);
	virtual _PMathObj CharAccess		(_PMathObj,_PMathObj);
	virtual _PMathObj Execute 			(long opCode, _MathObject* p = nil , _MathObject* p2 = nil);  
	virtual	BaseRef	  toStr 			(void);
	
	virtual	bool	  IsVariable		(void)
					  					{ return true; }
					  					
	virtual	bool	  HasChanged		(void)
					  					{ return true; }
	
	_String*	 	  theString;
		
};

//__________________________________________________________________________________

	
extern		_List 			variablePtrs;

extern		_SimpleList		BuiltInFunctionParameterCount,
							*deferSetFormula;

extern		_AVLListX		variableNames;

extern		_String			UnOps,
							HalfOps;

extern		_Parameter  	printDigits,
							verbosityLevel;

extern		long 	   		lastMatrixDeclared,
							subNumericValues;

long 		LocateVarByName (_String&);
_Variable*  LocateVar 		(long index);
_Variable*  FetchVar 		(long index);
_PMathObj	FetchObjectFromVariableByType 
							(_String*, int);
_PMathObj	FetchObjectFromVariableByTypeIndex 
							(long, int);
_String&    AppendContainerName
							(_String&, _VariableContainer*);

void 		DeleteVariable 	(_String&, bool deleteself = true);
void 		DeleteTreeVariable 	
							(_String&, _SimpleList&,bool);
void 		checkParameter 	(_String& name, _Parameter& dest, _Parameter def, _VariableContainer* = nil);
void 		stashParameter 	(_String& name, _Parameter  newVal, bool);
void 		setParameter   	(_String& name, _Parameter def, _String* = nil);
void 		setParameter   	(_String& name, _PMathObj  def, bool = true, _String* = nil);

long	 	VerbosityLevel (void);
void		ReplaceVar 		(_Variable*);
void		RenameVariable 	(_String*,_String*);
void		CompileListOfUserExpressions 
							(_SimpleList&,_List&, bool doAll = false);

void  		FindUnusedObjectName 
							(_String&, _String&, _List&, bool = false);

void  		FindUnusedObjectName 
							(_String&, _String&, _AVLListX&, bool = false);
							
bool		ExpressionCalculator 
							(void);
							
_Variable*  CheckReceptacle
							(_String*,_String, bool = true, bool = false);
							
bool	    CheckReceptacleAndStore
							(_String*,_String, bool, _PMathObj, bool = true);

void		FinishDeferredSF(void);
							
void		SetupOperationLists (void);
void		ExportIndVariables
							(_String&, _String&, _SimpleList*);
void		ExportDepVariables
							(_String&, _String&, _SimpleList*);
void		ExportCatVariables
							(_String&, _SimpleList*);
							
void		SplitVariablesIntoClasses
							(_SimpleList&, _SimpleList&, _SimpleList&, _SimpleList&);
bool		CheckEqual		(_Parameter,_Parameter);

extern		_AVLListX		_hyApplicationGlobals;



_Parameter	AddNumbers  (_Parameter, _Parameter);
_Parameter	SubNumbers  (_Parameter, _Parameter);
_Parameter	MultNumbers (_Parameter, _Parameter);
_Parameter	DivNumbers  (_Parameter, _Parameter);
_Parameter  EqualNumbers(_Parameter, _Parameter);
_Parameter  LessThan	(_Parameter, _Parameter);
_Parameter  GreaterThan (_Parameter, _Parameter);
_Parameter  LessThanE	(_Parameter, _Parameter);
_Parameter  GreaterThanE(_Parameter, _Parameter);
_Parameter	Power	    (_Parameter, _Parameter);
_Parameter	RandomNumber(_Parameter, _Parameter);
_Parameter	ExpNumbers  (_Parameter);
_Parameter	LogNumbers  (_Parameter);
_Parameter	MinusNumber (_Parameter);
_Parameter	MaxNumbers  (_Parameter, _Parameter);
_Parameter	MinNumbers  (_Parameter, _Parameter);
_Parameter  FastMxAccess(Ptr, _Parameter);

void		PopulateArraysForASimpleFormula
						(_SimpleList&, _SimpleFormulaDatum*);

void		WarnNotDefined (_PMathObj, long);

extern		_Parameter  pi_const;
extern		bool		useGlobalUpdateFlag;

#endif
