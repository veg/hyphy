#ifndef     __STACK__
#define     __STACK__

#include "mathobj.h"

//__________________________________________________________________________________
class   _Stack   //computational stack
{

    friend class _Formula;
    friend class _Operation;
    friend class _Variable;

public:

    _Stack (void);
    ~_Stack (void);

    bool      Push (_PMathObj);     // push object onto the stack
    _PMathObj Pop (bool del = true);            // pop object from the top of the stack
    long      StackDepth (void);    // returns the depth of the stack
    void      Reset (void);         // clear the stack

    virtual   void    Initialize (void);
    virtual   void    Duplicate (BaseRef);

protected:

    _List  theStack;
};

#endif
