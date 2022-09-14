/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-2006
 Primary Development:
 Sergei L Kosakovsky Pond (sergeilkp@mac.com)
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdfrost@ucsd.edu)
 Art FY Poon    (apoon@biomail.ucsd.edu)
 
 Some of the Original Code by William A Casey.
 
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

  // iterator implementation

template <class node_data> node_iterator<node_data>::node_iterator(
                                                                   node<node_data>* root, int traverser_kind){
  
  
  this->travseral_kind = traverser_kind;
  this->traversal_level      = -1L; // not initialized yet
  this->iterator_state       = root;
}


template <class node_data> long node_iterator<node_data>::Level(void) const {
  return this->iterator_state ? this->traversal_level : -1L;
}

template <class node_data> void node_iterator<node_data>::Reset(node<node_data>* root) {
  this->traversal_level      = -1L; // not initialized yet
  this->iterator_state       = root;
}

template <class node_data> node<node_data>* node_iterator<node_data>::Current(void) const {
  if (this->traversal_level >= 0L) {
    return this->iterator_state;
  }
  return NULL;
}

template <class node_data> void node_iterator<node_data>::pop_history_item (_SimpleList* history) {
  if (history)
    history->Pop();
}

template <class node_data> void node_iterator<node_data>::push_history_item (_SimpleList* history) {
  if (history)
    (*history)<< (long)this->iterator_state;
}

template <class node_data> node<node_data>* node_iterator<node_data>::Next(_SimpleList* history) {
  
  node<node_data>* test_node;
  
  if (this->iterator_state) {
    switch (this->travseral_kind) {
      case _HY_TREE_TRAVERSAL_POSTORDER:
      case _HY_TREE_TRAVERSAL_POSTORDER_RIGHT_FIRST: {
        
        if (this->traversal_level < 0L) {
          this->traversal_level = 0L;
          this->push_history_item (history);
          while ((test_node = this->iterator_state->go_down(this->travseral_kind == _HY_TREE_TRAVERSAL_POSTORDER ? 1 : iterator_state->get_num_nodes()))) {
            this->iterator_state = test_node;
            this->push_history_item (history);
            this->traversal_level++;
          }
          //printf ("[node_iterator init]%ld\n", this->traversal_level);
          return  this->iterator_state;
        }
        
        if (this->traversal_level == 0L) {
          // need this if traversal is initiated at an interior node which is not root
          this->iterator_state = nil;
        } else {
          test_node = this->iterator_state;
          test_node = this->travseral_kind == _HY_TREE_TRAVERSAL_POSTORDER ? test_node->go_next() : test_node->go_previous();
         
          if (test_node) {
            this->pop_history_item (history);
            this->iterator_state=test_node;
            this->push_history_item (history);
            while ((test_node = this->iterator_state->go_down(this->travseral_kind == _HY_TREE_TRAVERSAL_POSTORDER ? 1 : iterator_state->get_num_nodes()))) {
              this->iterator_state = test_node;
              this->push_history_item (history);
              this->traversal_level++;
            }
            return this->iterator_state;
          }
          this->iterator_state=this->iterator_state->go_up();
          this->pop_history_item (history);
          //printf ("[node_iterator up]%ld\n", this->traversal_level);
          if (this->traversal_level == 0L) {
            // need this if traversal is initiated at an interior node which is not root
            this->iterator_state = nil;
          } else {
            this->traversal_level--;
          }
        }
      }
      break;
        
      case _HY_TREE_TRAVERSAL_PREORDER: {
        
        if (this->traversal_level < 0L) {
          this->traversal_level = 0L;
        } else {
          
          test_node = this->iterator_state->go_down (1);
            // iterator_state has already been visited;
            // try to descend down the tree first
            // failing that, move to the sibling
            // failing that, move up the tree
          
          if (test_node) {
            this->iterator_state = test_node;
            this->traversal_level  ++;
            
          } else {
            test_node = this->iterator_state->go_next();
            if (test_node) {
              this->pop_history_item (history);
              this->iterator_state = test_node;
            } else {
              test_node = this->iterator_state->go_up();
              this->pop_history_item (history);
              while (test_node && test_node->go_down (test_node->get_num_nodes()) == this->iterator_state) {
                this->iterator_state = test_node;
                this->pop_history_item (history);
                this->traversal_level --;             
                test_node = this->iterator_state->go_up();
              }
              if (test_node) {
                this->iterator_state = this->iterator_state->go_next();
              } else {
                this->iterator_state = nil;
                break;
              }
            }
          }
        }
        this->push_history_item (history);
     }
      break;
    }
  }

  return this->iterator_state;
}

  //-------------------------------------------------------------
template <class node_data> void node<node_data>::delete_tree(bool delSelf){
 	
    long 	nc = get_num_nodes();
    for (int i=1; i<=nc; i++) {
        go_down(i)->delete_tree();
        delete (go_down(i));
    }
    if (delSelf) {
        delete (this);
    }
}

  //-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::duplicate_tree(void (callback) (node<node_data>*, node<node_data>*)) {
    
    node<node_data>* result = new node<node_data>;
    
    for (int i=1; i<=get_num_nodes(); i++) {
        result->add_node(*(go_down(i)->duplicate_tree(callback)));
    }
    

    if (callback) {
        callback (this, result);
    } else {
        result->in_object = in_object;
    }
    
    return result;
}

  //-------------------------------------------------------------

template <class node_data> void node_count_descendants (node<node_data>* , node<node_data>* n) {
    if (n->get_num_nodes() == 0) {
        n->in_object = 1L;
    } else {
        n->in_object = 0L;
        for (int i=1; i<=n->get_num_nodes(); i++) {
          n->in_object += n->go_down(i)->in_object;
        }
    }
}


  //-------------------------------------------------------------
template <class node_data> bool node<node_data>::compare_subtree(node<node_data>* compareTo)
{
  int nNodes = get_num_nodes();
  if (nNodes==compareTo->get_num_nodes()) {
    for (int i=1; i<=nNodes; i++)
      if (!go_down(i)->compare_subtree (compareTo->go_down(i)))
        return false;
    return true;
  }
  return false;
}


  //-------------------------------------------------------------

/* template <class node_data> long NodePathTraverser (_SimpleList& history, node<node_data>* root)
{
  static long  going_up, branchCount, tipCount;
  static node<node_data>* laststep;
  node<node_data>* curstep, *crashdummy;
		
  if (root)
  {
    laststep = root;
    branchCount=-1;
    tipCount=-1;
    history.Clear();
    while ((crashdummy = laststep->go_down(1)))
    {
      laststep = crashdummy;
      if (branchCount>-1)
        history<<branchCount;
      branchCount++;
    }
    tipCount=0;
    branchCount--;
    return 0;
  }
  
  curstep = laststep;
  crashdummy = curstep->go_next();
  if (crashdummy)
  {
    curstep=crashdummy;
    while ((crashdummy = curstep->go_down(1)))
    {
      branchCount++;
      history<<branchCount;
      curstep = crashdummy;
    }
    laststep = curstep;
    going_up = false;
    laststep = curstep;
    return ++tipCount;
  }
  
  curstep = curstep->parent;
  history.Delete(history.countitems()-1);
  if (!curstep) return -1;
  crashdummy = curstep->go_next();
  while (!crashdummy)
  {
    curstep=curstep->parent;
    if (!curstep)
    {
      if (!crashdummy) return -1;
    }
    crashdummy = curstep->go_next();
    history.Delete(history.countitems()-1);
  }
  going_up = true;
  laststep = curstep;
  return NodePathTraverser (history,(node<node_data>*)nil);
  
  return -1;
} */



  //-------------------------------------------------------------
template <class node_data> int node<node_data>::tree_depth(void)
{
	 int res = 0;
	 for (int i=nodes.get_length(); i>0;i--)
   {
     int t = go_down(i)->tree_depth();
     if (t>res)
       res = t;
   }
	 return res+1;
}

  //-------------------------------------------------------------
template <class node_data>
node<node_data>* NodeTraverser  (node<node_data>* root)
{
  static node<node_data>* laststep;
  node<node_data>* curstep, *crashdummy;
		
  if (root)
  {
    laststep = root;
    while ((crashdummy = laststep->go_down(1))) laststep = crashdummy;
    return laststep;
  }
  
  curstep = laststep;
  crashdummy = curstep->go_next();
  if (crashdummy)
  {
    curstep=crashdummy;
    while ((crashdummy = curstep->go_down(1))) curstep = crashdummy;
    return laststep = curstep;
  }
  curstep=curstep->get_parent();
  laststep = curstep;
  return curstep;
}
  //-----------------------------------Set Number 1----------------


template <class node_data> void node<node_data>::replace_node(node<node_data>* existing, node<node_data>* newNode){
    if (one == existing) {
        one = newNode;
        return;
    } else {
        if (two == existing) {
            two = newNode;
            return;
        }
    }
    for (long j = 0; nodes.length; j++) {
        if (nodes.data[j] == existing) {
          nodes.data[j] = newNode;
          break;
        }
    }
}

  //-----------------------------------Set Number 1----------------
template <class node_data> int node<node_data>::get_child_num() {
    if (parent) {
        if (this == parent->one) return 1;
        if (this == parent->two) return 2;
        for (int i=0; i<parent->nodes.length; i++) {
            if (parent->nodes.data[i] == this) return (i+3);
        }
    }
    return -1;
}
  //----------------------------------end no 1--------------------
  //---------Public set no 2-------------------------------------

  //--Bool (T/F) responses to potential tree moves
/*template <class node_data> int node<node_data>::down(int index)
{
  if ((index > 0) && (index <= get_num_nodes())){
    return 1;
  }
  else return 0;
} //Truth of decent*/

  //-------------------------------------------------------------

/*template <class node_data> int node<node_data>::next(){
  
  if (get_parent() == NULL) return 0;
  if (get_child_num() < (get_parent())->get_num_nodes()) return 1;
  return 0;
}

  //-------------------------------------------------------------
template <class node_data> int node<node_data>::up(){
  if (get_child_num() > 0) return 1;
  else return 0;
}*/

  //--------MOVERS MOVERS through the tree
template <class node_data> node<node_data>* node<node_data>::go_up(){
  return get_parent();
}

  //-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_next(){
    
    int my_index = get_child_num();
    if (my_index > 0) {
        if (my_index < parent->get_num_nodes()) {
            return parent->get_node (my_index + 1);
        }
    }
    return NULL;
}

  //-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_previous(){
    int my_index = get_child_num();
    if (my_index > 1) {
        return parent->get_node (my_index - 1);
    }
    return NULL;
}


  //-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_down(int index){
  if (index > 0L && index <= get_num_nodes())
  	 return (get_node(index));
  
  return NULL; //false=can go no further
}
  //-------------end public set no 2-----------------------------
