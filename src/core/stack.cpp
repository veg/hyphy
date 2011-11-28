#include "stack.h"

//__________________________________________________________________________________
_Stack::_Stack (void)
{
}

//__________________________________________________________________________________
void _Stack::Initialize (void)
{
    theStack.Initialize();
}

//__________________________________________________________________________________
void _Stack::Duplicate (BaseRef s)
{
    theStack.Duplicate(&((_Stack*)s)->theStack);
}

//__________________________________________________________________________________
_Stack::~_Stack (void)
{
}

//__________________________________________________________________________________
bool _Stack::Push (_PMathObj newObj)    // push object onto the stack
{
    theStack<<(newObj);
    return true;
}

//__________________________________________________________________________________
_PMathObj _Stack::Pop (bool del)        // pop object from the top of the stack
{
    _PMathObj r = (_PMathObj)theStack.lData[theStack.lLength-1];
    if (del) {
        theStack.lLength--;
    }
    return r;
}

//__________________________________________________________________________________
long _Stack::StackDepth (void)  // returns the depth of the stack
{
    return theStack.lLength;
}

//__________________________________________________________________________________
void _Stack::Reset (void)   // clears the stack
{
    theStack.Clear();
}
