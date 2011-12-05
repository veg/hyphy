#ifndef     __OPERATION__
#define     __OPERATION__

#include "baseobj.h"
#include "list.h"
#include "hy_strings.h"
#include "mathobj.h"

extern  _List BuiltInFunctions;


class _Stack;
class _VariableContainer;
class _Formula;
//__________________________________________________________________________________
class   _Operation : public BaseObj
{

    friend class _Formula;
    friend class _Variable;
    friend class _VariableContainer;

public:
    _Operation  (void);
    _Operation  (_String&, long);
    // construct the operation by its symbol and, if relevant -
    // number of operands

    _Operation  (bool, _String&, bool isG = false, _VariableContainer*  = nil);
    // store a variable or a constant
    _Operation  (_PMathObj);
    // store a non-numeric constant

    virtual ~_Operation (void);

    virtual   BaseObj*      makeDynamic         (void);

    bool            Execute             (_Stack&, _VariableContainer* = nil); //execute this operation
    // see the commend for _Formula::ExecuteFormula for the second argument
    virtual   void          StackDepth          (long&);

    bool            ExecutePolynomial   (_Stack&);
    virtual   BaseObj*      toStr               (void);    //convert the op to string

    virtual   void          Initialize          (void);
    virtual   void          Duplicate           (BaseRef);
    _String&    GetCode             (void) {
        return (opCode>-1)&&(numberOfTerms>=0)?*(_String*)BuiltInFunctions(opCode):empty;
    }
    long&       TheCode             (void) {
        return opCode;
    }
    virtual  bool           IsAVariable         (bool = true) ; // is this object a variable or not?
    virtual  bool           IsConstant          (void);         // does this object depend on any independent variables or not?

    virtual  long           UserFunctionID      (void) {
        return numberOfTerms < 0 ? -numberOfTerms-1 : -1;
    };
    // return a non-neg number (function index) if this is a user function,
    // otherwise, return -1

    virtual  long           GetAVariable        (void) {    // return the index of the variable
        return theData>=-2?theData:-theData-3;
    }

    virtual  void           SetAVariable        (long d) {  // return the index of the variable
        theData=d;
    }

    virtual  bool           AssignmentVariable  (void) {
        return theData<-2;
    }

    virtual  bool           HasChanged          (void);

    virtual  void           SetTerms            (long d) {
        numberOfTerms=d;
    }

    virtual  _PMathObj      GetANumber          (void) {
        return theNumber;
    }

    virtual  void           SetNumber           (_PMathObj d) {
        theNumber=d;
    }

    long            GetNoTerms          (void) {
        return numberOfTerms;
    }
    long            PrecedenceLevel     (void);

    bool            CanResultsBeCached (_Operation *);


    virtual bool            EqualOp             (_Operation*);

protected:


    long        opCode;         // internal operation code
    long        numberOfTerms,  // 1 - unary, 2 - binary, etc
                theData;
    _PMathObj   theNumber;
};

#endif
