/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
  Steven Weaver (sweaver@ucsd.edu)
  
Module Developers:
	Lance Hepler (nlhepler@gmail.com)
	Martin Smith (martin.audacis@gmail.com)

Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdf22@cam.ac.uk)

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


//#pragma once
#ifndef __CLASSES__
#define __CLASSES__

#ifndef MAX
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#define SWAP(a,b,c) ((c)=(a),(a)=(b),(b)=(c))
#define EXCHANGE(a,b) {auto c = (a); (a) = (b); (b) = c;} 
#define DEBUG 0

#include "simplelist.h"

/*---------------------------------------------------------------------------*/

template <class array_data> class ptr_array {
    template<typename T> friend class node;
protected:
    int             length;                                            //length of the array
    
public:
    array_data  *data;                                            //points to an array of somethings
 
    ptr_array                   (void) : data (NULL), length (0) {
    }
    //default constructor
    ~ptr_array              () {
        if (data) {
            delete [] data;
        }
    }
    //default destructor
    void            add                 (array_data);
    //adds a node to the end of the array
    void            prepend             (array_data);
    //adds a node to the beginning of the array
    void            delete_entry        (int);
    void            clear               () {
      if (data) {
        delete [] data;
        data = NULL;
        length = 0;
      }
    };
    int             get_length          (void) const {
        return length;
    }
    //returns the length of the array  array_data get_elem(int index);
};

/*---------------------------------------------------------------------------*/

template <class node_data> class node {
public:
    node_data       in_object;
    ptr_array<node*>    nodes;
    node                *parent,  *one,  *two;
    void                set_parent      (node &thenode) {
        parent = &thenode;
    }
public:
    node            (void) {
        parent = NULL;
        one = NULL;
        two = NULL;
    }
    node            (node_data const & d) {
      parent = NULL; one = NULL; two = NULL;
      this->init (d);
    }
    ~node()
    {}
    void                init            (const node_data & value) {
        in_object = value;
    }
    node_data       get_data        (void) {
        return in_object;
    }
    //TEST BELOW
    void                detach_child    (int k) {
      int do_delete = 0;
      if (k == 1 || k == 2) {
        if (k == 1) {
          one = two;
        }
        if (nodes.length) {
          two = nodes.data[0];
          do_delete = 1;
        } else {
          two = NULL;
        }
      } else {
        do_delete = k - 2;
      }
      if (do_delete > 0) {
          nodes.delete_entry(do_delete);
      }
    }

    int                 tree_depth      (void);

    void                detach_parent   (void) {
        parent = NULL;
    }
    
    int                 get_num_nodes   (void) const {
        if (one) {
            if (two) {
                return 2+nodes.get_length();
            }
            return 1;
        }
        return 0;
    }
    
    void                add_node        (node& thenode) {
        thenode.set_parent(*this);
        if (one) {
            if (two) {
                nodes.add(&thenode);
            } else {
                two = &thenode;
            }
        } else {
            one = &thenode;
        }
    }
    
    template <typename... Args> void add_node (node& first, Args&... args) {
        add_node (first);
        add_node (args...);
    }


    void                replace_node    (node* existingNode, node* newNode);

    void                prepend_node    (node &thenode) {
        thenode.set_parent(*this);
        if (!one) {
            one = &thenode;
        } else {
            if (two) {
                nodes.prepend (two);
            }
            two = one;
            one = &thenode;
        }
        //nodes.prepend(&thenode);
    }
    void                kill_node       (int in) {
        
      if (in < 3) {
          if (in == 1) {
              one = two;
          }
          if (nodes.length) {
              two = nodes.data[0];
              nodes.delete_entry(1);
          } else {
              two = NULL;
          }
      } else {
          nodes.delete_entry(in - 2);
      }
    }

  
    void                kill_all_nodes       (void) {
       one = NULL; two = NULL;
       nodes.clear();
    }
  
    node*           get_node        (int in) {
        if (in == 1) return one;
        if (in == 2) return two;
        return nodes.data[in-3];
    }
    node*           get_parent      (void) const {
        return parent;
    }
    int                 get_child_num   (void);
    bool                is_leaf (void) const {return get_num_nodes() == 0L; }
    bool                is_root (void) const {return get_parent() == nil;}
  
    node<node_data>*    go_up           (void);
    node<node_data>*    go_next           (void);
    node<node_data>*    go_previous       (void);
    node<node_data>*    go_down           (int index);
    int                 down            (int index);

    int                 next            (void);
    int                 up              (void);
    void                delete_tree     (bool = false);
    node<node_data>*    duplicate_tree  (void (callback) (node<node_data>*, node<node_data>*) = nil);
    // traverse from this node down duplicating the tree and return
    // pointer to the root of the duplicate;
    bool                compare_subtree (node<node_data>*);

};

#define _HY_TREE_TRAVERSAL_POSTORDER              0x00
#define _HY_TREE_TRAVERSAL_POSTORDER_RIGHT_FIRST  0x01
#define _HY_TREE_TRAVERSAL_PREORDER               0x02

template <class node_data> class node_iterator {
  private:
    node<node_data>*       iterator_state;
    long                   traversal_level;
    int                    travseral_kind;
    void                   push_history_item  (_SimpleList* history);
    void                   pop_history_item   (_SimpleList* history);
  
  public:
    node_iterator   (node<node_data>* root, int traverser_kind = _HY_TREE_TRAVERSAL_POSTORDER);
    ~node_iterator  () {};
    void             Reset (node<node_data>* root);
    node<node_data>* Next (_SimpleList* history = NULL);
    node<node_data>* Current (void) const;
    long             Level (void) const;
};


#include "../ptr_array.cp"
#include "../classes.cp"

#endif
