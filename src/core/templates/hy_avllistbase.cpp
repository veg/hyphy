/*

 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Steven Weaver (sweaver@ucsd.edu)
 Martin Smith (martin.audacis@gmail.com)
 
 Module Developers:
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Lance Hepler (nlhepler@gmail.com)
 
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

    _AVLList    structure inspired by the excellent documentation of
    GNU libavl 2.0.1 by Ben Pfaff (http://www.msu.edu/~pfaffben/avl/index.html)

*/

//*************** CONSTRUCTORS ***************//


template <typename KEYTYPE, typename PAYLOAD>
_AVLListBase<KEYTYPE,PAYLOAD>::_AVLListBase(void) {
    root = HY_AVL_LEAF;
}

template <typename KEYTYPE, typename PAYLOAD>
_AVLListBase<KEYTYPE,PAYLOAD>::_AVLListBase(_AVLListBase<KEYTYPE,PAYLOAD> const & source) {
    //keys.Clone (&source.keys);
    //values.Clone (&source.values);
    avl_structure.Clone (&source.avl_structure);
    root = source.root;
}


//*************** DESTRUCTOR ***************//

template <typename KEYTYPE, typename PAYLOAD>
_AVLListBase<KEYTYPE,PAYLOAD>::~_AVLListBase(void) {
    
}


//*************** ATTRIBUTE ACCESSORS ***************//

template <typename KEYTYPE, typename PAYLOAD>
unsigned long _AVLListBase<KEYTYPE,PAYLOAD>::Length (void) const {
    return this->avl_structure.Length() - this->empty_slots.Length();
}


//*************** SEARCH FUNCTIONS  ***************//

template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Find (const KEYTYPE& key, _AVLListTraversalHistory * history) const{
    long current_node = this->root;
    
    while (current_node != HY_AVL_LEAF) {
        long comp = this->_CompareIndexToValue (current_node, key);
        if (comp) {
            if (history) {
                history->push (current_node);
            }
            current_node = this->_MoveInTree (current_node, comp);
         } else {
            return current_node;
        }
    }
    
    return HY_NOT_FOUND;
}

template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::FindBest (const KEYTYPE& key, long &last_index ) const{
    long current_node = root,
         comp = HY_AVL_MOVE_RIGHT;
    
    while (current_node != HY_AVL_LEAF) {
        comp = this->_CompareIndexToValue (current_node, key);
        last_index = current_node;
        
        if (comp == 0L) {
            return 0L;
        }

        current_node = this->_MoveInTree (current_node, comp);

    }
    
    return comp;
}


//*************** HELPER FUNCTIONS ***************//

template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::_MoveInTree (long node, long direction) const{
    if (direction < 0L) {
        return this->avl_structure.AtIndex(node).left_child;
    } else {
        if (direction > 0L) {
            return this->avl_structure.AtIndex(node).right_child;
        }
    }
        
    return node;
}


//*************** TREE TRAVERSAL FUNCTIONS ***************//

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::_DescendToTerminal(long current_index, long direction, _AVLListTraversalHistory * history) const {
    
    long the_index = HY_NOT_FOUND;
    
    while (current_index != HY_AVL_LEAF) {
        if (history) {
            history->push (current_index);
        }
        the_index = current_index;
        current_index = this->_MoveInTree (current_index, direction);
    }
    
    return the_index;
}


//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::First(_AVLListTraversalHistory * history) const {
    return this->_DescendToTerminal (this->root, HY_AVL_MOVE_LEFT, history);
}

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Last(_AVLListTraversalHistory * history) const {
    return this->_FindExtremeValue (this->root, HY_AVL_MOVE_RIGHT, history);
}

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Next(_AVLListTraversalHistory& history) const {
    
    if (history.Length()) {
        long current_node = history.pop(),
             try_node = this->_MoveInTree (current_node, HY_AVL_MOVE_RIGHT);
        if (try_node != HY_AVL_LEAF) {
            return this->_DescentToTerminal (try_node, HY_AVL_MOVE_LEFT, &history);
        } else {
             while (history.Length()) {
                long parent_node = history.pop();
                try_node = this->_MoveInTree (parent_node, HY_AVL_MOVE_RIGHT);
                if (try_node == current_node) {
                    current_node = parent_node;
                } else {
                    history.push (try_node);
                    return try_node;
                }
            }
        }
    }
    
    return HY_NOT_FOUND;
}

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Prev(_AVLListTraversalHistory& history) const {
    
    if (history.Length()) {
        long current_node = history.pop(),
        try_node = this->_MoveInTree (current_node, HY_AVL_MOVE_LEFT);
        if (try_node != HY_AVL_LEAF) {
            return this->_DescentToTerminal (try_node, HY_AVL_MOVE_RIGHT, &history);
        } else {
            while (history.Length()) {
                long parent_node = history.pop();
                try_node = this->_MoveInTree (parent_node, HY_AVL_MOVE_LEFT);
                if (try_node == current_node) {
                    current_node = parent_node;
                } else {
                    history.push (try_node);
                    return try_node;
                }
            }
        }
    }
    
    return HY_NOT_FOUND;
}

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Traverser(_AVLListTraversalHistory &history) const {
    if (history.Length()) {
        return this->Next (&history);
    } else {
        return this->First (&history);
    }
    return HY_NOT_FOUND;
}


//*************** INSERTION AND DELETION FUNCTIONS ***************//

//______________________________________________________________________________
template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::Insert(KEYTYPE const& key, PAYLOAD const& value) {
    if (this->Length() > 0UL) { // something to do
        
        
        _SimpleList directions;
        
        long y = root,
             z = HY_AVL_LEAF,
             p,
             q,
             n,
             w;
        
        //bool go_right = false;
        int move_direction = HY_AVL_MOVE_LEFT;
        
        for (q = z, p = y; p != HY_AVL_LEAF;
             q = p, p = this->_MoveInTree (p, move_direction)) {
            
            //long comp = dataList->Compare(b, p);
            long move_direction = this->_CompareIndexToValue (p, key);
            
            if (move_direction == 0L) {
                 return -p - 1L;
            }
 
            if (this->avl_structure[p].balance_factor != 0L) {
                z = q;
                y = p;
                directions.Clear();
            }
            directions.append (move_direction);
        }
        
        
        // insert new node
        
        n = this->_StoreKeyValuePair (key, value);
        
        if (move_direction == HY_AVL_MOVE_RIGHT) {
            this->avl_structure[q].right_child = n;
        } else {
            this->avl_structure[q].left_child  = n;
        }
        
        // update balance factors
        
        p = y;
        
        for (long k = 0L; p != n;
             p =  this->_MoveInTree (p, directions.AtIndex(k)), k++) {
            this->avl_structure.AtIndex(k).balance_factor += directions.AtIndex(k) == HY_AVL_MOVE_RIGHT ? 1L : -1L;
        
        }
        
        if (avl_structure.AtIndex(y).balance_factor  == -2L) {
            long x = this->avl_structure.AtIndex(y).left_child;
            if (this->avl_structure.AtIndex(x).balance_factor == -1) {
                w = x;
                LEFT_SHIFT (this->avl_structure.AtIndex(y).left_child,this->avl_structure.AtIndex(x).right_child,y);
                this->avl_structure.AtIndex(x).balance_factor = 0L;
                this->avl_structure.AtIndex(y).balance_factor = 0L;
            } else {
                w = this->avl_structure.AtIndex(x).right_child;
                LEFT_SHIFT (this->avl_structure.AtIndex(x).right_child,this->avl_structure.AtIndex(w).left_child,x);
                LEFT_SHIFT (this->avl_structure.AtIndex(y).left_child,this->avl_structure.AtIndex(w).right_child,y);
                
                switch (this->avl_structure.AtIndex(w).balance_factor) {
                    case -1L:
                        this->avl_structure.AtIndex(x).balance_factor = 0L;
                        this->avl_structure.AtIndex(y).balance_factor = 1L;
                        this->avl_structure.AtIndex(w).balance_factor = 0L;
                        break;
                    case 0L:
                        this->avl_structure.AtIndex(x).balance_factor = 0L;
                        this->avl_structure.AtIndex(y).balance_factor = 0L;
                        break;
                    default:
                        this->avl_structure.AtIndex(x).balance_factor = -1L;
                        this->avl_structure.AtIndex(y).balance_factor = 0L;
                        this->avl_structure.AtIndex(w).balance_factor = 0L;
                        break;
                }
                
            }
        } else if (avl_structure.AtIndex(y).balance_factor == 2L) {
            long x = this->avl_structure.AtIndex(y).right_child;
            if (this->avl_structure.AtIndex(x).balance_factor == 1L) {
                w = x;
                LEFT_SHIFT (this->avl_structure.AtIndex(y).right_child, this->avl_structure.AtIndex(x).left_child,y);
                this->avl_structure.AtIndex(x).balance_factor = 0L;
                this->avl_structure.AtIndex(y).balance_factor = 0L;
            } else {
                w = this->avl_structure.AtIndex(x).left_child;
                LEFT_SHIFT (this->avl_structure.AtIndex(x).left_child,this->avl_structure.AtIndex(w).right_child,x);
                LEFT_SHIFT (this->avl_structure.AtIndex(y).right_child,this->avl_structure.AtIndex(w).left_child,y);
                switch (this->avl_structure.AtIndex(w).balance_factor) {
                    case 1L:
                        this->avl_structure.AtIndex(x).balance_factor = 0L;
                        this->avl_structure.AtIndex(y).balance_factor = -1L;
                        this->avl_structure.AtIndex(w).balance_factor = 0L;
                        break;
                    case 0L:
                        this->avl_structure.AtIndex(x).balance_factor = 0L;
                        this->avl_structure.AtIndex(y).balance_factor = 0L;
                        break;
                    default:
                        this->avl_structure.AtIndex(x).balance_factor = 1L;
                        this->avl_structure.AtIndex(y).balance_factor = 0L;
                        this->avl_structure.AtIndex(w).balance_factor = 0L;
                        break;
                }
               
            }
        } else {
            return n;
        }
        
        if (z != HY_AVL_LEAF) {
            if (y == this->avl_structure.AtIndex(z).left_child) {
                this->avl_structure.AtIndex(z).left_child = w;
            } else {
                this->avl_structure.AtIndex(z).right_child = w;
            }
        }
        
        if (y == this->root) {
            this->root = w;
        }        
        return p;
    }
    
    
    return (this->root = this->_StoreKeyValuePair (HY_LIST_INSERT_AT_END, key, value));
    
}

/*template <typename KEYTYPE,typename PAYLOAD>
long _AVLListBase<KEYTYPE,PAYLOAD>::_CompareIndexToValue (long node, KEYTYPE const & key) const{
    return this->keys.CompareToValue (node, key);
}*/

/*






//______________________________________________________________________________
long _AVLList::GetByIndex(const long theIndex) {
  if (theIndex == 0) {
    return First();
  }

  long elementCount = countitems();

  if (theIndex == elementCount - 1) {
    return Last();
  }

  if (theIndex > 0 && theIndex < elementCount) {
    _SimpleList hist;
    long ls, cn = Traverser(hist, ls, GetRoot()), counter = 0;

    while (counter < theIndex) {
      counter++;
      cn = Traverser(hist, ls);
    }

    return cn;

  }

  return HY_NOT_FOUND;
}


//______________________________________________________________________________
void _AVLList::ReorderList(_SimpleList *s) {
  _SimpleList reorderMe(
      (unsigned long)(dataList->lLength - emptySlots.lLength + 1)),
      nodeStack((unsigned long) 32);

  long curNode = root;

  while (1) {
    while (curNode >= 0) {
      nodeStack << curNode;
      curNode = leftChild.lData[curNode];
    }
    if (long h = nodeStack.lLength) {
      h--;
      curNode = nodeStack.lData[h];
      if (s) {
        (*s) << curNode;
      }
      reorderMe.InsertElement(((BaseRef *)dataList->lData)[curNode], -1, false,
                              false);

      curNode = rightChild.lData[curNode];
      nodeStack.Delete(h, false);

    } else {
      break;
    }
  }

  reorderMe.TrimMemory();

  long *t = dataList->lData;
  dataList->lData = reorderMe.lData;
  dataList->lLength = reorderMe.lLength;
  dataList->laLength = reorderMe.laLength;
  reorderMe.lData = t;
}

//______________________________________________________________________________
void _AVLList::ConsistencyCheck(void) {
  _SimpleList nodeStack((unsigned long) 32);

  long curNode = root, lastNode = -1, checkCount = 0;

  while (1) {
    while (curNode >= 0) {
      nodeStack << curNode;
      curNode = leftChild.lData[curNode];
      if (curNode >= (long) dataList->lLength) {
        WarnError("Failed Constistency Check in _AVLList");
        return;
      }

    }
    if (long h = nodeStack.lLength) {
      if (h > 3 * log(1. + countitems())) {
        WarnError("Failed Constistency Check in _AVLList");
        return;
      }
      h--;
      curNode = nodeStack.lData[h];
      if (lastNode >= 0 && curNode >= 0) {
        if (dataList->Compare(Retrieve(lastNode), curNode) >= 0) {
          WarnError("Failed Constistency Check in _AVLList");
          return;
        }
        checkCount++;
      }
      if ((balanceFactor.lData[curNode] < -1) ||
          (balanceFactor.lData[curNode] > 1)) {
        WarnError("Failed Constistency Check in _AVLList");
        return;
      }
      lastNode = curNode;
      curNode = rightChild.lData[curNode];
      if (curNode >= (long) dataList->lLength) {
        WarnError("Failed Constistency Check in _AVLList");
        return;
      }
      nodeStack.Delete(h, false);
    } else {
      break;
    }
  }

  if (dataList->lLength &&
      (dataList->lLength > checkCount + 1 + emptySlots.lLength)) {
    WarnError("Failed Constistency Check in _AVLList");
    return;
  }

}


//______________________________________________________________________________
BaseRef _AVLList::toStr(void) {
  _String *str = new _String(128L, true);
  checkPointer(str);

  if (countitems() == 0) {
    (*str) << "Empty Associative List";
  } else {
    _SimpleList hist;
    long ls, cn;

    cn = Traverser(hist, ls, root);

    while (cn >= 0) {
      long keyVal = (long) Retrieve(cn);
      (*str) << _String(keyVal);
      (*str) << '\n';
      cn = Traverser(hist, ls);
    }
  }

  str->Finalize();
  return str;
}

//______________________________________________________________________________
BaseRef _AVLList::Retrieve(long idx) {
  return ((BaseRef *)dataList->lData)[idx];
}

//______________________________________________________________________________
void _AVLList::Clear(bool cL) {
  if (cL) {
    dynamic_cast<_List*>(dataList)->Clear();
  } else {
    dataList->Clear();
  }

  emptySlots.Clear();
  root = -1;
  leftChild.Clear();
  rightChild.Clear();
  balanceFactor.Clear();
}

//______________________________________________________________________________
long _AVLList::InsertData(BaseRef b, long, bool) {
  long w = (long) emptySlots.lLength - 1, n;

  if (w >= 0) {
    n = emptySlots.lData[w];
    emptySlots.Delete(w);
    leftChild.lData[n] = -1;
    rightChild.lData[n] = -1;
    balanceFactor.lData[n] = 0;
    ((BaseRef *)dataList->lData)[n] = b;
  } else {
    n = dataList->lLength;
    dataList->InsertElement(b, -1, false, false);
    leftChild << -1;
    rightChild << -1;
    balanceFactor << 0;
  }
  return n;
}


//______________________________________________________________________________
long _AVLList::Insert(BaseRef b, long xtra, bool cp, bool clear) {
  if (dataList->lLength - emptySlots.lLength) {
    long y = root, z = -1, p, q, n, w;

    bool go_right = false;

    _SimpleList da((unsigned long) 32);

    for (q = z, p = y; p >= 0;
         q = p, p = go_right ? rightChild.lData[p] : leftChild.lData[p]) {
      long comp = dataList->Compare(b, p);
      if (comp == 0) {
        if (cp == false && clear) {
          DeleteObject(b);
        }
        return -p - 1;
      }
      if (balanceFactor.lData[p] != 0) {
        z = q;
        y = p;
        da.Clear();
      }
      go_right = comp > 0;
      da << go_right;
    }


    // insert new node

    n = InsertData(b, xtra, cp);

    if (go_right) {
      rightChild.lData[q] = n;
    } else {
      leftChild.lData[q] = n;
    }

    // update balance factors

    p = y;

    for (long k = 0; p != n;
         p = da.lData[k] ? rightChild.lData[p] : leftChild.lData[p], k++)
      if (da.lData[k] == 0) {
        balanceFactor.lData[p]--;
      } else {
        balanceFactor.lData[p]++;
      }

    //if (z < 0)
    //{
    //ConsistencyCheck();
    //return n;
    //}

    if (balanceFactor.lData[y] == -2) {
      //152

      long x = leftChild.lData[y];
      if (balanceFactor.lData[x] == -1) { //155
        w = x;
        leftChild.lData[y] = rightChild.lData[x];
        rightChild.lData[x] = y;
        balanceFactor.lData[x] = balanceFactor.lData[y] = 0;
      } else { //156
        w = rightChild.lData[x];
        rightChild.lData[x] = leftChild.lData[w];
        leftChild.lData[w] = x;
        leftChild.lData[y] = rightChild.lData[w];
        rightChild.lData[w] = y;
        if (balanceFactor.lData[w] == -1) {
          balanceFactor.lData[x] = 0;
          balanceFactor.lData[y] = 1;
        } else if (balanceFactor.lData[w] == 0) {
          balanceFactor.lData[x] = 0;
          balanceFactor.lData[y] = 0;
        } else {
          balanceFactor.lData[x] = -1;
          balanceFactor.lData[y] = 0;
        }

        balanceFactor.lData[w] = 0;
      }
    } else if (balanceFactor.lData[y] == 2) {
      long x = rightChild.lData[y];
      if (balanceFactor.lData[x] == 1) {
        w = x;
        rightChild.lData[y] = leftChild.lData[x];
        leftChild.lData[x] = y;
        balanceFactor.lData[x] = balanceFactor.lData[y] = 0;
      } else {
        w = leftChild.lData[x];
        leftChild.lData[x] = rightChild.lData[w];
        rightChild.lData[w] = x;
        rightChild.lData[y] = leftChild.lData[w];
        leftChild.lData[w] = y;
        if (balanceFactor.lData[w] == 1) {
          balanceFactor.lData[x] = 0;
          balanceFactor.lData[y] = -1;
        } else if (balanceFactor.lData[w] == 0) {
          balanceFactor.lData[x] = 0;
          balanceFactor.lData[y] = 0;
        } else {
          balanceFactor.lData[x] = 1;
          balanceFactor.lData[y] = 0;
        }

        balanceFactor.lData[w] = 0;
      }
    } else {
      //ConsistencyCheck ();
      return n;
    }

    if (z >= 0) {
      if (y == leftChild.lData[z]) {
        leftChild.lData[z] = w;
      } else {
        rightChild.lData[z] = w;
      }
    }

    if (y == root) {
      root = w;
    }

    //ConsistencyCheck ();

    return p;
  }


  root = InsertData(b, xtra, cp);

  return 0;

}

//______________________________________________________________________________
bool _AVLList::HasData(long idx) { return leftChild.lData[idx] != 2; }

//______________________________________________________________________________
void _AVLList::Delete(BaseRef b, bool delMe) {

  if (root == -1) {
    return;
  }

  _SimpleList pa((unsigned long) 64), da((unsigned long) 64);

  long p = root, cmp = dataList->Compare(b, p), k = 0;

  pa.lData[k] = -1;
  da.lData[k++] = 1;

  for (; cmp != 0; cmp = dataList->Compare(b, p)) {
    bool go_right = cmp > 0;

    pa.lData[k] = p;
    da.lData[k++] = go_right;

    if (go_right) {
      p = rightChild.lData[p];
    } else {
      p = leftChild.lData[p];
    }

    if (p < 0) {
      return;
    }
  }

  if (k == 1) {
    pa.lData[k] = -1;
  }

  emptySlots << p;
  if (delMe) {
    DeleteObject(Retrieve(p));
  }
  //((BaseRef*)dataList->lData)[p] = nil;
  dataList->lData[p] = 0;
  DeleteXtra(p);

  long r = rightChild.lData[p];

  if (r < 0) {
    if (k > 1) {
      if (da.lData[k - 1] == 1) {
        rightChild.lData[pa.lData[k - 1]] = leftChild.lData[p];
      } else {
        leftChild.lData[pa.lData[k - 1]] = leftChild.lData[p];
      }
    }

    if (p == root)
        //root = pa.lData[k-1];
        {
      root = leftChild.lData[root];
    }
  } else {
    if (leftChild.lData[r] < 0) {
      leftChild.lData[r] = leftChild.lData[p];
      balanceFactor.lData[r] = balanceFactor.lData[p];
      if (k > 1) {
        if (da.lData[k - 1] == 1) {
          rightChild.lData[pa.lData[k - 1]] = r;
        } else {
          leftChild.lData[pa.lData[k - 1]] = r;
        }
      } else {
        root = r;
      }

      da.lData[k] = 1;
      pa.lData[k++] = r;
      //if (p==root)
      //root = r;
    } else {
      long s;
      int j = k++;
      for (;;) {
        da.lData[k] = 0;
        pa.lData[k++] = r;
        s = leftChild.lData[r];
        if (leftChild.lData[s] < 0) {
          break;
        }
        r = s;
      }

      leftChild.lData[s] = leftChild.lData[p];
      leftChild.lData[r] = rightChild.lData[s];
      rightChild.lData[s] = rightChild.lData[p];
      balanceFactor.lData[s] = balanceFactor.lData[p];

      if (j > 1) {
        if (da.lData[j - 1] == 1) {
          rightChild.lData[pa.lData[j - 1]] = s;
        } else {
          leftChild.lData[pa.lData[j - 1]] = s;
        }
      }

      da.lData[j] = 1;
      pa.lData[j] = s;
      if (p == root) {
        root = s;
      }
    }
  }

  //if (k>63)
  //{
  //WarnError ("Internal List error");
  //}

  while (--k > 0) {
    long y = pa.lData[k];
    if (da.lData[k] == 0) {
      balanceFactor.lData[y]++;
      if (balanceFactor.lData[y] == 1) {
        break;
      } else if (balanceFactor.lData[y] == 2) {
        long x = rightChild.lData[y];
        if (balanceFactor.lData[x] == -1) {
          long w = leftChild.lData[x];
          leftChild.lData[x] = rightChild.lData[w];
          rightChild.lData[w] = x;
          rightChild.lData[y] = leftChild.lData[w];
          leftChild.lData[w] = y;
          if (balanceFactor.lData[w] == 1) {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = -1;
          } else if (balanceFactor.lData[w] == 0) {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = 0;
          } else {
            balanceFactor.lData[x] = 1;
            balanceFactor.lData[y] = 0;
          }

          balanceFactor.lData[w] = 0;
          if (k > 1) {
            if (da.lData[k - 1] == 1) {
              rightChild.lData[pa.lData[k - 1]] = w;
            } else {
              leftChild.lData[pa.lData[k - 1]] = w;
            }
          } else {
            root = w;
          }

          //if (y==root)
          //  root = w;
        } else {
          rightChild.lData[y] = leftChild.lData[x];
          leftChild.lData[x] = y;

          if (k > 1) {
            if (da.lData[k - 1] == 1) {
              rightChild.lData[pa.lData[k - 1]] = x;
            } else {
              leftChild.lData[pa.lData[k - 1]] = x;
            }
          } else {
            root = x;
          }

          if (balanceFactor.lData[x] == 0) {
            balanceFactor.lData[x] = -1;
            balanceFactor.lData[y] = 1;
            break;
          } else {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = 0;
          }
        }
      }
    } else {
      balanceFactor.lData[y]--;
      if (balanceFactor.lData[y] == -1) {
        break;
      } else if (balanceFactor.lData[y] == -2) {
        long x = leftChild.lData[y];
        if (balanceFactor.lData[x] == 1) {
          long w = rightChild.lData[x];
          rightChild.lData[x] = leftChild.lData[w];
          leftChild.lData[w] = x;
          leftChild.lData[y] = rightChild.lData[w];
          rightChild.lData[w] = y;
          if (balanceFactor.lData[w] == -1) {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = 1;
          } else if (balanceFactor.lData[w] == 0) {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = 0;
          } else {
            balanceFactor.lData[x] = -1;
            balanceFactor.lData[y] = 0;
          }

          balanceFactor.lData[w] = 0;
          if (k > 1) {
            if (da.lData[k - 1] == 1) {
              rightChild.lData[pa.lData[k - 1]] = w;
            } else {
              leftChild.lData[pa.lData[k - 1]] = w;
            }
          } else {
            root = w;
          }

          //if (y==root)
          //root = w;
        } else {
          leftChild.lData[y] = rightChild.lData[x];
          rightChild.lData[x] = y;
          if (k > 1) {
            if (da.lData[k - 1] == 1) {
              rightChild.lData[pa.lData[k - 1]] = x;
            } else {
              leftChild.lData[pa.lData[k - 1]] = x;
            }
          } else {
            root = x;
          }

          if (balanceFactor.lData[x] == 0) {
            balanceFactor.lData[x] = 1;
            balanceFactor.lData[y] = -1;
            break;
          } else {
            balanceFactor.lData[x] = 0;
            balanceFactor.lData[y] = 0;
          }
        }
      }
    }
  }
  //ConsistencyCheck ();

} */
