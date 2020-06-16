/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 */


#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include "avllist.h"
#include "hy_strings.h"
#include "hy_string_buffer.h"
#include "parser.h"

#include "global_things.h"

using namespace hy_global;


//______________________________________________________________
// AVL Lists
//______________________________________________________________

_AVLList::_AVLList (_SimpleList* d) {
    dataList = d;
    root     = -1;
}

_AVLList::_AVLList (_AVLList const &) {
    HandleApplicationError("Called _AVLList:: copy constructor  method stub that is not implemented");
}

//______________________________________________________________

BaseRef _AVLList::makeDynamic (void) const {
    HandleApplicationError("Called _AVLList::makeDynamic:  method stub that is not implemented");
    return nil;
}

//______________________________________________________________

void _AVLList::Duplicate (BaseRefConst) {
    HandleApplicationError("Called _AVLList::Duplicate:  method stub that is not implemented");
}

//______________________________________________________________

void _AVLList::operator = (_AVLList const & rhs) {
    HandleApplicationError("Called _AVLList::operator = :  method stub that is not implemented");
}

//______________________________________________________________

long  _AVLList::Find (BaseRefConst obj) const {
    long curNode = root;

    while (curNode>=0) {
        long comp = dataList->Compare (obj,curNode);

        if (comp<0) {
            curNode = leftChild.list_data[curNode];
        } else if (comp>0) {
            curNode = rightChild.list_data[curNode];
        } else {
            return curNode;
        }
    }

    return kNotFound;
}

//______________________________________________________________

long  _AVLList::FindLong (long obj) const {
    long curNode = root;

    while (curNode>=0) {
        long comp = dataList->list_data[curNode];

        if (obj<comp) {
            curNode = leftChild.list_data[curNode];
        } else if (obj>comp) {
            curNode = rightChild.list_data[curNode];
        } else {
            return curNode;
        }
    }

    return kNotFound;
}

//______________________________________________________________

char  _AVLList::FindBest (BaseRefConst obj, long& lastNode) const {
    long curNode  = root,
         comp     = 1;

    while (curNode>=0 && comp) {
        comp = dataList->Compare (obj,curNode);
        lastNode = curNode;

        if (comp<0) {
            curNode = leftChild.list_data[curNode];
        } else if (comp>0) {
            curNode = rightChild.list_data[curNode];
        } else {
            return 0;
        }
    }

    return comp;
}

//______________________________________________________________

long  _AVLList::Find (BaseRefConst obj, _SimpleList& hist) const {
    long curNode = root;

    while (curNode>=0) {
        long comp = dataList->Compare (obj,curNode);

        if (comp<0) {
            hist << curNode;
            curNode = leftChild.list_data[curNode];
        } else if (comp>0) {
            hist << curNode;
            curNode = rightChild.list_data[curNode];
        } else {
            return curNode;
        }
    }

    return -1;
}

//______________________________________________________________

long  _AVLList::Next (long d, _SimpleList& hist) const {
    if (d >= 0) {
        if (rightChild.list_data [d] >= 0) {
            hist << d;
            d = rightChild.list_data [d];
            while (leftChild.list_data[d] >= 0) {
                hist << d;
                d = leftChild.list_data[d];
            }
            return d;
        } else {
            while (hist.countitems()) {
                long x = hist.Pop();

                if (rightChild.list_data[x] != d) {
                    return x;
                }
                //TODO:???
                d = x;
            }

            return -1;
        }
    }

    d = root;
    while (d >= 0 && leftChild.list_data[d] >=0) {
      hist << d;
      d = leftChild.list_data[d];
    }
    
    return d;
}

//______________________________________________________________

long  _AVLList::First (void) const {
    long   d = root;
    while (d >= 0 && leftChild.list_data[d] >=0) {
        d = leftChild.list_data[d];
    }

    return d;
}

//______________________________________________________________

long  _AVLList::Last (void) const {
    long   d = root;
    while (d >= 0 && rightChild.list_data[d] >=0) {
        d = rightChild.list_data[d];
    }

    return d;
}

//______________________________________________________________

long  _AVLList::Prev (long d, _SimpleList& hist) const {
  if (d >= 0) {
    if (leftChild.list_data [d] >= 0) {
      hist << d;
      d = leftChild.list_data [d];
      while (rightChild.list_data[d] >= 0) {
        hist << d;
        d = rightChild.list_data[d];
      }
      return d;
    } else {
      while (hist.countitems()) {
        long x = hist.list_data[hist.lLength-1];
        
        hist.Delete (hist.lLength-1);
        
        if (leftChild.list_data[x] != d) {
          return x;
        }
        //TODO:???
        d = x;
      }
      
      return -1;
    }
  }
  
  d = root;
  while (d >= 0 && rightChild.list_data[d] >=0) {
    hist << d;
    d = rightChild.list_data[d];
  }
  
  return d;
  
}


//______________________________________________________________

bool  _AVLList::IsValidIndex(long index) const {
  if (index >= 0 && index < dataList->lLength) {
    return Retrieve(index);
  }
  return false;
}

//______________________________________________________________

long  _AVLList::GetByIndex (const long theIndex) {
    if (theIndex == 0) {
        return First();
    }

    long elementCount = countitems();

    if (theIndex == elementCount - 1) {
        return Last();
    }

    if (theIndex > 0 && theIndex < elementCount) {
        _SimpleList  hist;
        long         ls,
                     cn = Traverser (hist,ls,GetRoot()),
                     counter = 0;

        while (counter < theIndex) {
            counter ++;
            cn = Traverser (hist,ls);
        }

        return cn;

    }

    return -1;
}

//______________________________________________________________

void  _AVLList::ReorderList (_SimpleList *s) {
    
    unsigned long item_count = (dataList->lLength-emptySlots.lLength+1);
    
    _SimpleList reorderMe (item_count),
                nodeStack (64UL, (long*)alloca (64*sizeof (unsigned long)));

    long        curNode = root;

    while (1) {
        while (curNode >= 0) {
            nodeStack << curNode;
            curNode = leftChild.list_data[curNode];
        }
        if (long h = nodeStack.lLength) {
            h--;
            curNode = nodeStack.list_data[h];
            if (s) {
                (*s) << curNode;
            }
            
            reorderMe.InsertElement (((BaseRef*)dataList->list_data)[curNode],-1,false,false);
            curNode = rightChild.list_data[curNode];
            nodeStack.Delete (h, false);

        } else {
            break;
        }
    }

    reorderMe.TrimMemory ();
    *dataList = reorderMe;
}

//______________________________________________________________

void  _AVLList::ConsistencyCheck (void)
{
    _SimpleList nodeStack ((unsigned long)32);

    long        curNode  = root,
                lastNode = -1,
                checkCount = 0;

    while (1) {
        while (curNode >= 0) {
            nodeStack << curNode;
            curNode = leftChild.list_data[curNode];
            if (curNode >= (long)dataList->lLength) {
                hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
                return;
            }

        }
        if (long h = nodeStack.lLength) {
            if (h>3*log (1.+countitems())) {
                hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
                return;
            }
            h--;
            curNode = nodeStack.list_data[h];
            if (lastNode >= 0 && curNode >= 0) {
                if (dataList->Compare (Retrieve (lastNode), curNode) >= 0) {
                    hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
                    return;
                }
                checkCount++;
            }
            if ((balanceFactor.list_data[curNode] < -1)||(balanceFactor.list_data[curNode] > 1)) {
                hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
                return;
            }
            lastNode = curNode;
            curNode = rightChild.list_data[curNode];
            if (curNode >= (long)dataList->lLength) {
                hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
                return;
            }
            nodeStack.Delete (h, false);
        } else {
            break;
        }
    }

    if (dataList->lLength && (dataList->lLength > checkCount + 1 + emptySlots.lLength)) {
        hy_global::HandleApplicationError ("Failed Constistency Check in _AVLList");
        return;
    }

}

//______________________________________________________________

long  _AVLList::Traverser (_SimpleList &nodeStack, long& t, long r) const {
    if  (r >= 0) {
        t = r;
        nodeStack.Clear();
    }

    while (t >= 0) {
        nodeStack << t;
        t = leftChild.list_data[t];
    }

    if (long h = nodeStack.lLength) {
        h--;
        t = nodeStack.list_data[h];
        r = t;
        t = rightChild.list_data[t];
        nodeStack.Delete (h, false);
        return r;
    }
    return -1;
}

//______________________________________________________________

BaseRef  _AVLList::toStr (unsigned long) {
    _StringBuffer * str = new _StringBuffer (128L);
 
    if (countitems() == 0) {
        (*str) << "()";
    } else {
        _SimpleList  hist;
        long         ls, cn;

        cn = Traverser (hist,ls,root);
      
        bool first = true;
      
        (*str) << '(';
      
        while (cn>=0) {
            if (first) {
              first = false;
            } else {
              (*str) << ", ";
            }
           (*str) << _String((long)Retrieve (cn));
            cn = Traverser (hist,ls);
        }
      
        (*str) << ')';
    }

    return str;
}

  //______________________________________________________________

const _List  _AVLList::Keys (void) const {
  _List keys;
   if (countitems() > 0UL) {
      _SimpleList  hist;
      long         ls = -1L, cn;
      
      cn = Traverser (hist,ls,root);
      
      while (cn>=0) {
        keys << Retrieve (cn);
        cn = Traverser (hist,ls);
      }
  }
  return keys;
}


//______________________________________________________________

BaseRef _AVLList::Retrieve (long idx) const {
    return ((BaseRef*)dataList->list_data)[idx];
}

  //______________________________________________________________

long _AVLList::RetrieveLong(long idx) const {
  return ((long*)dataList->list_data)[idx];
}

//______________________________________________________________

void _AVLList::Clear (bool cL)
{
    if (cL) {
        ((_List*)dataList)->Clear();
    } else {
        dataList->Clear();
    }

    emptySlots.Clear();
    root = -1;
    leftChild.Clear();
    rightChild.Clear();
    balanceFactor.Clear();
}

//______________________________________________________________

long  _AVLList::InsertData (BaseRef b, long, bool)
{
    long w = (long)emptySlots.lLength - 1,
         n;

    if (w>=0) {
        n = emptySlots.list_data[w];
        emptySlots.Delete (w);
        leftChild.list_data[n] = -1;
        rightChild.list_data[n] = -1;
        balanceFactor.list_data[n] = 0;
        ((BaseRef*)dataList->list_data)[n] = b;
    } else {
        n = dataList->lLength;
        dataList->InsertElement (b,-1,false,false);
        leftChild  << -1;
        rightChild << -1;
        balanceFactor << 0;
    }
    return n;
}

//______________________________________________________________

unsigned long _AVLList::countitems (void) const {
    return dataList->lLength - emptySlots.lLength;
}

//______________________________________________________________

long  _AVLList::Insert (BaseRef b, long xtra,bool cp,bool clear) {
/** 
  Insert a key (and possibly a value) into this _AVLList
 
  @param b the key to insert; this operation will NOT increase reference counts for b
  @param xtra the payload
  @param cp if true, the insertion operation will add a reference count to payload
  @param clear if insertion fails (key already exists), then the _key_ will be deleted
 
 */
    if (dataList->lLength-emptySlots.lLength) {
        long        y = root,
                    z = -1,
                    p,
                    q,
                    n,
                    w;

        bool        go_right = false;

        //_SimpleList da ((unsigned long)32);
        
        long          traversal_stack [32] = {0L};
        unsigned long stack_pointer = 0UL;

        // try to find the node or where to insert it

        for (q=z, p=y; p>=0; q=p, p=go_right?rightChild.list_data[p]:leftChild.list_data[p]) {
            long comp = dataList->Compare (b, p);
            if (comp == 0) {
                if (cp == false && clear) {
                    DeleteObject (b);
                }
                // already exists, return the node it mapped to
                return -p-1;
            }
            if (balanceFactor.list_data[p] != 0) {
                z = q;
                y = p;
                //da.Clear();
                stack_pointer = 0UL;
            }
            go_right = comp > 0;
            //da << go_right;
            traversal_stack[stack_pointer++] = go_right;
            if (stack_pointer >= 32UL) {
                HandleApplicationError("Internal error: AVLList is too big");
            }
        }

        /*if (da.lLength > 3*log (dataList->lLength+2))
        {
            WarnError ("AVLList internal error!");
            return -1;
        }*/

        // insert new node

        n = InsertData (b, xtra,cp);

        if (go_right) {
            rightChild.list_data[q] = n;
        } else {
            leftChild.list_data[q] = n;
        }


        // update balance factors

        p = y;

        for (long k=0L; p!=n; p=traversal_stack[k]?rightChild.list_data[p]:leftChild.list_data[p],k++)
            //if (da.list_data[k] == 0) {
            if (traversal_stack[k] == 0) {
                balanceFactor.list_data[p]--;
            } else {
                balanceFactor.list_data[p]++;
            }


        //if (z < 0)
        //{
        //ConsistencyCheck();
        //return n;
        //}

        if (balanceFactor.list_data[y] == -2) {
            //152

            long x = leftChild.list_data[y];
            if (balanceFactor.list_data[x] == -1) { //155
                w                    = x;
                leftChild.list_data [y]  = rightChild.list_data[x];
                rightChild.list_data[x]  = y;
                balanceFactor.list_data[x] = balanceFactor.list_data[y] = 0;
            } else { //156
                w = rightChild.list_data[x];
                rightChild.list_data[x] = leftChild.list_data[w];
                leftChild.list_data[w] = x;
                leftChild.list_data[y] = rightChild.list_data[w];
                rightChild.list_data[w] = y;
                if (balanceFactor.list_data[w] == -1) {
                    balanceFactor.list_data[x] = 0;
                    balanceFactor.list_data[y] = 1;
                } else if (balanceFactor.list_data[w] == 0) {
                    balanceFactor.list_data[x] = 0;
                    balanceFactor.list_data[y] = 0;
                } else {
                    balanceFactor.list_data[x] = -1;
                    balanceFactor.list_data[y] = 0;
                }

                balanceFactor.list_data[w] = 0;
            }
        } else if (balanceFactor.list_data[y] == 2) {
            long x = rightChild.list_data[y];
            if (balanceFactor.list_data[x] == 1) {
                w                      = x;
                rightChild.list_data [y]   = leftChild.list_data[x];
                leftChild.list_data[x]     = y;
                balanceFactor.list_data[x] = balanceFactor.list_data[y] = 0;
            } else {
                w = leftChild.list_data[x];
                leftChild.list_data[x] = rightChild.list_data[w];
                rightChild.list_data[w] = x;
                rightChild.list_data[y] = leftChild.list_data[w];
                leftChild.list_data[w] = y;
                if (balanceFactor.list_data[w] == 1) {
                    balanceFactor.list_data[x] = 0;
                    balanceFactor.list_data[y] = -1;
                } else if (balanceFactor.list_data[w] == 0) {
                    balanceFactor.list_data[x] = 0;
                    balanceFactor.list_data[y] = 0;
                } else {
                    balanceFactor.list_data[x] = 1;
                    balanceFactor.list_data[y] = 0;
                }

                balanceFactor.list_data[w] = 0;
            }
        } else {
            //ConsistencyCheck ();
            return n;
        }

        if (z >= 0) {
            if (y == leftChild.list_data[z]) {
                leftChild.list_data[z] = w;
            } else {
                rightChild.list_data[z] = w;
            }
        }

        if (y==root) {
            root = w;
        }

        //ConsistencyCheck ();

        return p;
    }

    root =InsertData (b, xtra,cp);

    return 0;

}

//______________________________________________________________

bool  _AVLList::HasData (long idx) {
    return leftChild.list_data[idx] != 2;
}

//______________________________________________________________

void  _AVLList::Delete (BaseRefConst b, bool delMe) {

    if (root == -1) {
        return;
    }

    _SimpleList pa ((unsigned long)64),
                da ((unsigned long)64);

    long        p = root,
                cmp = dataList->Compare (b,p),
                k = 0;

    pa.list_data[k] = -1;
    da.list_data[k++] = 1;

    for (; cmp !=0; cmp = dataList->Compare (b,p)) {
        bool go_right = cmp > 0;

        pa.list_data[k] = p;
        da.list_data[k++] = go_right;

        if (go_right) {
            p = rightChild.list_data[p];
        } else {
            p = leftChild.list_data[p];
        }

        if (p<0) {
            return;
        }
    }

    if (k==1) {
        pa.list_data[k]   = -1;
    }

    emptySlots << p;
    if (delMe) {
        DeleteObject (Retrieve(p));
    }
    //((BaseRef*)dataList->list_data)[p] = nil;
    dataList->list_data[p] = 0;
    DeleteXtra (p);

    long r = rightChild.list_data[p];

    if (r < 0) {
        if (k>1) {
            if (da.list_data[k-1] == 1) {
                rightChild.list_data[pa.list_data[k-1]] = leftChild.list_data[p];
            } else {
                leftChild.list_data[pa.list_data[k-1]] = leftChild.list_data[p];
            }
        }

        if (p==root) {
            //root = pa.list_data[k-1];
            root = leftChild.list_data[root];
        }
    } else {
        if (leftChild.list_data[r] < 0) {
            leftChild.list_data[r]     = leftChild.list_data[p];
            balanceFactor.list_data[r] = balanceFactor.list_data[p];
            if (k>1) {
                if (da.list_data[k-1] == 1) {
                    rightChild.list_data[pa.list_data[k-1]] = r;
                } else {
                    leftChild.list_data[pa.list_data[k-1]] = r;
                }
            } else {
                root = r;
            }

            da.list_data[k]   = 1;
            pa.list_data[k++] = r;
            //if (p==root)
            //root = r;
        } else {
            long s;
            int  j = k++;
            for (;;) {
                da.list_data[k]   = 0;
                pa.list_data[k++] = r;
                s = leftChild.list_data[r];
                if (leftChild.list_data[s] < 0) {
                    break;
                }
                r = s;
            }


            leftChild.list_data[s] = leftChild.list_data[p];
            leftChild.list_data[r] = rightChild.list_data[s];
            rightChild.list_data[s] = rightChild.list_data[p];
            balanceFactor.list_data[s] = balanceFactor.list_data[p];

            if (j>1) {
                if (da.list_data[j-1] == 1) {
                    rightChild.list_data[pa.list_data[j-1]] = s;
                } else {
                    leftChild.list_data[pa.list_data[j-1]] = s;
                }
            }

            da.list_data[j] = 1;
            pa.list_data[j] = s;
            if (p==root) {
                root = s;
            }
        }
    }

    //if (k>63)
    //{
    //WarnError ("Internal List error");
    //}

    while (--k > 0) {
        long y = pa.list_data[k];
        if (da.list_data[k] == 0) {
            balanceFactor.list_data[y] ++;
            if (balanceFactor.list_data[y] == 1) {
                break;
            } else if (balanceFactor.list_data[y] == 2) {
                long x = rightChild.list_data[y];
                if (balanceFactor.list_data[x] == -1) {
                    long w = leftChild.list_data[x];
                    leftChild.list_data[x] = rightChild.list_data[w];
                    rightChild.list_data[w] = x;
                    rightChild.list_data[y] = leftChild.list_data[w];
                    leftChild.list_data[w] = y;
                    if (balanceFactor.list_data[w] == 1) {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = -1;
                    } else if (balanceFactor.list_data[w] == 0) {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = 0;
                    } else {
                        balanceFactor.list_data[x] = 1;
                        balanceFactor.list_data[y] = 0;
                    }

                    balanceFactor.list_data[w] = 0;
                    if (k>1) {
                        if (da.list_data[k-1] == 1) {
                            rightChild.list_data[pa.list_data[k-1]] = w;
                        } else {
                            leftChild.list_data[pa.list_data[k-1]] = w;
                        }
                    } else {
                        root = w;
                    }

                    //if (y==root)
                    //  root = w;
                } else {
                    rightChild.list_data[y] = leftChild.list_data[x];
                    leftChild.list_data[x] = y;

                    if (k>1) {
                        if (da.list_data[k-1] == 1) {
                            rightChild.list_data[pa.list_data[k-1]] = x;
                        } else {
                            leftChild.list_data[pa.list_data[k-1]] = x;
                        }
                    } else {
                        root = x;
                    }

                    if (balanceFactor.list_data[x] == 0) {
                        balanceFactor.list_data[x] = -1;
                        balanceFactor.list_data[y] = 1;
                        break;
                    } else {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = 0;
                    }
                }
            }
        } else {
            balanceFactor.list_data[y] --;
            if (balanceFactor.list_data[y] == -1) {
                break;
            } else if ( balanceFactor.list_data[y] == -2) {
                long x = leftChild.list_data[y];
                if (balanceFactor.list_data[x] == 1) {
                    long w = rightChild.list_data[x];
                    rightChild.list_data[x] = leftChild.list_data[w];
                    leftChild.list_data[w] = x;
                    leftChild.list_data[y] = rightChild.list_data[w];
                    rightChild.list_data[w] = y;
                    if (balanceFactor.list_data[w] == -1) {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = 1;
                    } else if (balanceFactor.list_data[w] == 0) {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = 0;
                    } else {
                        balanceFactor.list_data[x] = -1;
                        balanceFactor.list_data[y] = 0;
                    }

                    balanceFactor.list_data[w] = 0;
                    if (k>1) {
                        if (da.list_data[k-1] == 1) {
                            rightChild.list_data[pa.list_data[k-1]] = w;
                        } else {
                            leftChild.list_data[pa.list_data[k-1]] = w;
                        }
                    } else {
                        root = w;
                    }

                    //if (y==root)
                    //root = w;
                } else {
                    leftChild.list_data[y] = rightChild.list_data[x];
                    rightChild.list_data[x] = y;
                    if (k>1) {
                        if (da.list_data[k-1] == 1) {
                            rightChild.list_data[pa.list_data[k-1]] = x;
                        } else {
                            leftChild.list_data[pa.list_data[k-1]] = x;
                        }
                    } else {
                        root = x;
                    }

                    if (balanceFactor.list_data[x] == 0) {
                        balanceFactor.list_data[x] = 1;
                        balanceFactor.list_data[y] = -1;
                        break;
                    } else {
                        balanceFactor.list_data[x] = 0;
                        balanceFactor.list_data[y] = 0;
                    }
                }
            }
        }
    }
    //ConsistencyCheck ();

}
