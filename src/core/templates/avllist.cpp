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
_AVLList<KEYTYPE,PAYLOAD>::_AVLList(void) {
    root = HY_AVL_LEAF;
}

template <typename KEYTYPE, typename PAYLOAD>
_AVLList<KEYTYPE,PAYLOAD>::_AVLList(_AVLList<KEYTYPE,PAYLOAD> const & source) {
    keys.Clone (&source.keys);
    values.Clone (&source.values);
    avl_structure.Clone (&source.avl_structure);
    root = source.root;
}


//*************** DESTRUCTOR ***************//

template <typename KEYTYPE, typename PAYLOAD>
_AVLList<KEYTYPE,PAYLOAD>::~_AVLList(void) {
    
}


//*************** ATTRIBUTE ACCESSORS ***************//

template <typename KEYTYPE, typename PAYLOAD>
unsigned long _AVLList<KEYTYPE,PAYLOAD>::Length (void) const {
    return this->avl_structure.Length() - this->empty_slots.Length();
}


//*************** SEARCH FUNCTIONS  ***************//

template <typename KEYTYPE,typename PAYLOAD>
long _AVLList<KEYTYPE,PAYLOAD>::Find (const KEYTYPE key, _SimpleList* history) const{
    long current_node = this->root;
    
    while (current_node != HY_AVL_LEAF) {
        long comp = this->_CompareIndexToValue (current_node, key);
        if (comp) {
            if (history) {
                history->append(current_node);
            }
            current_node = this->_MoveInTree (current_node, comp);
         } else {
            return current_node;
        }
    }
    
    return HY_NOT_FOUND;
}

//*************** HELPER FUNCTIONS ***************//

template <typename KEYTYPE,typename PAYLOAD>
long _AVLList<KEYTYPE,PAYLOAD>::_MoveInTree (long node, long direction) const{
    if (direction < 0L) {
        return this->avl_structure(node).left_child;
    } else {
        if (direction > 0L) {
            return this->avl_structure(node).right_child;
        }
    }
        
    return node;
}

template <typename KEYTYPE,typename PAYLOAD>
long _AVLList<KEYTYPE,PAYLOAD>::_CompareIndexToValue (long node, KEYTYPE const & key) const{
    return this->keys.CompareToValue (node, key);
}
/*


//______________________________________________________________________________
char _AVLList::FindBest(BaseRef obj, long &lastNode) {
  long curNode = root, comp = 1;

  while (curNode >= 0 && comp) {
    comp = dataList->Compare(obj, curNode);
    lastNode = curNode;

    if (comp < 0) {
      curNode = leftChild.lData[curNode];
    } else if (comp > 0) {
      curNode = rightChild.lData[curNode];
    } else {
      return 0;
    }
  }

  return comp;
}

//______________________________________________________________________________
long _AVLList::Find(BaseRef obj, _SimpleList &hist) {
  long curNode = root;

  while (curNode >= 0) {
    long comp = dataList->Compare(obj, curNode);

    if (comp < 0) {
      hist << curNode;
      curNode = leftChild.lData[curNode];
    } else if (comp > 0) {
      hist << curNode;
      curNode = rightChild.lData[curNode];
    } else {
      return curNode;
    }
  }

  return HY_NOT_FOUND;
}

//______________________________________________________________________________
long _AVLList::Next(long d, _SimpleList &hist) {
  if (d >= 0) {
    if (rightChild.lData[d] >= 0) {
      hist << d;
      d = rightChild.lData[d];
      while (leftChild.lData[d] >= 0) {
        hist << d;
        d = leftChild.lData[d];
      }
      return d;
    } else {
      while (hist.countitems()) {
        long x = hist.lData[hist.lLength - 1];

        hist.Delete(hist.lLength - 1);

        if (rightChild.lData[x] != d) {
          return x;
        }
        d = x;
      }

      return HY_NOT_FOUND;
    }
  }

  d = root;
  while (d >= 0 && leftChild.lData[d] >= 0) {
    d = leftChild.lData[d];
  }

  return d;
}

//______________________________________________________________________________
long _AVLList::First(void) {
  long d = root;
  while (d >= 0 && leftChild.lData[d] >= 0) {
    d = leftChild.lData[d];
  }

  return d;
}

//______________________________________________________________________________
long _AVLList::Last(void) {
  long d = root;
  while (d >= 0 && rightChild.lData[d] >= 0) {
    d = rightChild.lData[d];
  }

  return d;
}

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
long _AVLList::Prev(long d, _SimpleList &hist) {
  if (d >= 0) {
    if (leftChild.lData[d] >= 0) {
      hist << d;
      d = leftChild.lData[d];
      while (rightChild.lData[d] >= 0) {
        hist << d;
        d = rightChild.lData[d];
      }
      return d;
    } else {
      while (hist.countitems()) {
        long x = hist.lData[hist.lLength - 1];

        hist.Delete(hist.lLength - 1);

        if (leftChild.lData[x] != d) {
          return x;
        }
        d = x;
      }

      return HY_NOT_FOUND;
    }
  }

  d = root;
  while (d >= 0 && rightChild.lData[d] >= 0) {
    d = rightChild.lData[d];
  }

  return d;

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
long _AVLList::Traverser(_SimpleList &nodeStack, long &t, long r) {
  if (r >= 0) {
    t = r;
    nodeStack.Clear();
  }

  while (t >= 0) {
    nodeStack << t;
    t = leftChild.lData[t];
  }

  if (long h = nodeStack.lLength) {
    h--;
    t = nodeStack.lData[h];
    r = t;
    t = rightChild.lData[t];
    nodeStack.Delete(h, false);
    return r;
  }
  return HY_NOT_FOUND;
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
unsigned long _AVLList::countitems(void) {
  return dataList->lLength - emptySlots.lLength;
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
