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

// #pragma once
#ifndef __CLASSES__
#define __CLASSES__

#ifndef MAX
#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define SWAP(a, b, c) ((c) = (a), (a) = (b), (b) = (c))
#define EXCHANGE(a, b)                                                         \
  {                                                                            \
    auto c = (a);                                                              \
    (a) = (b);                                                                 \
    (b) = c;                                                                   \
  }
#define DEBUG 0

#include "simplelist.h"

/*---------------------------------------------------------------------------*/

/**
 * @brief A pointer array class
 *
 * @tparam array_data
 */
template <class array_data> class ptr_array {
  template <typename T> friend class node;

protected:
  int length; // length of the array

public:
  array_data *data; // points to an array of somethings

  /**
   * @brief Construct a new ptr_array object
   *
   */
  ptr_array(void) : length(0), data(NULL) {}
  // default constructor
  /**
   * @brief Destroy the ptr_array object
   *
   */
  ~ptr_array() {
    if (data) {
      delete[] data;
    }
  }
  // default destructor
  /**
   * @brief Add a node to the end of the array.
   *
   * @param ad The data to add.
   */
  void add(array_data);
  /**
   * @brief Add a node to the beginning of the array.
   *
   * @param ad The data to add.
   */
  void prepend(array_data);
  /**
   * @brief Delete an entry from the array.
   *
   * @param i The index of the entry to delete.
   */
  void delete_entry(int);
  /**
   * @brief Clear the array
   *
   */
  void clear() {
    if (data) {
      delete[] data;
      data = NULL;
      length = 0;
    }
  };
  /**
   * @brief Get the length of the array
   *
   * @return int
   */
  int get_length(void) const { return length; }
  // returns the length of the array  array_data get_elem(int index);
};

/*---------------------------------------------------------------------------*/

/**
 * @brief A node class
 *
 * @tparam node_data
 */
template <class node_data> class node {
public:
  node_data in_object;
  ptr_array<node *> nodes;
  node *parent, *one, *two;
  /**
   * @brief Set the parent of the node
   *
   * @param thenode
   */
  void set_parent(node &thenode) { parent = &thenode; }

public:
  /**
   * @brief Construct a new node object.
   */
  node(void) {
    parent = NULL;
    one = NULL;
    two = NULL;
  }
  /**
   * @brief Construct a new node object with data.
   *
   * @param d The data to initialize the node with.
   */
  node(node_data const &d) {
    parent = NULL;
    one = NULL;
    two = NULL;
    this->init(d);
  }
  /**
   * @brief Destroy the node object.
   */
  ~node() {}
  /**
   * @brief Initialize the node with data.
   *
   * @param value The data to initialize the node with.
   */
  void init(const node_data &value) { in_object = value; }
  /**
   * @brief Get the data of the node.
   *
   * @return The data of the node.
   */
  node_data get_data(void) { return in_object; }
  /**
   * @brief Detach a child from the node.
   *
   * @param k The index of the child to detach.
   */
  void detach_child(int k) {
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

  /**
   * @brief Get the depth of the tree starting from this node.
   *
   * @return The depth of the tree.
   */
  int tree_depth(void);

  /**
   * @brief Detach the parent of the node.
   */
  void detach_parent(void) { parent = NULL; }

  /**
   * @brief Get the number of children of this node.
   *
   * @return The number of children.
   */
  int get_num_nodes(void) const {
    if (one) {
      if (two) {
        return 2 + nodes.get_length();
      }
      return 1;
    }
    return 0;
  }

  /**
   * @brief Add a child node to this node.
   *
   * @param thenode The child node to add.
   */
  void add_node(node &thenode) {
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

  /**
   * @brief Add multiple child nodes to this node.
   *
   * @tparam Args The types of the other child nodes.
   * @param first The first child node to add.
   * @param args The other child nodes to add.
   */
  template <typename... Args> void add_node(node &first, Args &...args) {
    add_node(first);
    add_node(args...);
  }

  /**
   * @brief Replace a child node with another node.
   *
   * @param existingNode The child node to replace.
   * @param newNode The new node.
   */
  void replace_node(node *existingNode, node *newNode);

  /**
   * @brief Prepend a child node to this node.
   *
   * @param thenode The child node to prepend.
   */
  void prepend_node(node &thenode) {
    thenode.set_parent(*this);
    if (!one) {
      one = &thenode;
    } else {
      if (two) {
        nodes.prepend(two);
      }
      two = one;
      one = &thenode;
    }
    // nodes.prepend(&thenode);
  }
  /**
   * @brief Remove a child node from this node.
   *
   * @param in The index of the child node to remove.
   */
  void kill_node(int in) {

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

  /**
   * @brief Remove all child nodes from this node.
   */
  void kill_all_nodes(void) {
    one = NULL;
    two = NULL;
    nodes.clear();
  }

  /**
   * @brief Get a child node by its index.
   *
   * @param in The index of the child node.
   * @return A pointer to the child node.
   */
  node *get_node(int in) {
    if (in == 1)
      return one;
    if (in == 2)
      return two;
    return nodes.data[in - 3];
  }
  /**
   * @brief Get the parent of this node.
   *
   * @return A pointer to the parent node.
   */
  node *get_parent(void) const { return parent; }
  /**
   * @brief Get the index of this node in its parent's child list.
   *
   * @return The index of this node.
   */
  int get_child_num(void);
  /**
   * @brief Check if this node is a leaf node.
   *
   * @return true if this node is a leaf, false otherwise.
   */
  bool is_leaf(void) const { return get_num_nodes() == 0L; }
  /**
   * @brief Check if this node is the root of the tree.
   *
   * @return true if this node is the root, false otherwise.
   */
  bool is_root(void) const { return get_parent() == nil; }

  /**
   * @brief Go up to the parent node.
   *
   * @return A pointer to the parent node.
   */
  node<node_data> *go_up(void);
  /**
   * @brief Go to the next sibling node.
   *
   * @return A pointer to the next sibling node.
   */
  node<node_data> *go_next(void);
  /**
   * @brief Go to the previous sibling node.
   *
   * @return A pointer to the previous sibling node.
   */
  node<node_data> *go_previous(void);
  /**
   * @brief Go down to a child node.
   *
   * @param index The index of the child node.
   * @return A pointer to the child node.
   */
  node<node_data> *go_down(int index);
  /**
   * @brief Go down to a child node.
   *
   * @param index The index of the child node.
   * @return The index of the child node.
   */
  int down(int index);

  /**
   * @brief Go to the next sibling node.
   *
   * @return The index of the next sibling node.
   */
  int next(void);
  /**
   * @brief Go up to the parent node.
   *
   * @return The index of the parent node.
   */
  int up(void);
  /**
   * @brief Delete the subtree starting from this node.
   *
   * @param b If true, the node itself will be deleted.
   */
  void delete_tree(bool = false);
  /**
   * @brief Duplicate the subtree starting from this node.
   *
   * @param callback A callback function to be called for each duplicated node.
   * @return A pointer to the root of the duplicated subtree.
   */
  node<node_data> *duplicate_tree(void(callback)(node<node_data> *,
                                                 node<node_data> *) = nil);
  /**
   * @brief Compare this subtree with another subtree.
   *
   * @param n The subtree to compare with.
   * @return true if the subtrees are identical, false otherwise.
   */
  bool compare_subtree(node<node_data> *);
};

#define _HY_TREE_TRAVERSAL_POSTORDER 0x00
#define _HY_TREE_TRAVERSAL_POSTORDER_RIGHT_FIRST 0x01
#define _HY_TREE_TRAVERSAL_PREORDER 0x02

/**
 * @brief A node iterator class
 *
 * @tparam node_data
 */
template <class node_data> class node_iterator {
private:
  node<node_data> *iterator_state;
  long traversal_level;
  int travseral_kind;
  /**
   * @brief Push a history item to the history list.
   *
   * @param history The history list.
   */
  void push_history_item(_SimpleList *history);
  /**
   * @brief Pop a history item from the history list.
   *
   * @param history The history list.
   */
  void pop_history_item(_SimpleList *history);

public:
  /**
   * @brief Construct a new node_iterator object.
   *
   * @param root The root of the tree to iterate over.
   * @param traverser_kind The kind of traversal to perform.
   */
  node_iterator(node<node_data> *root,
                int traverser_kind = _HY_TREE_TRAVERSAL_POSTORDER);
  /**
   * @brief Destroy the node_iterator object.
   */
  ~node_iterator() {};
  /**
   * @brief Reset the iterator to the root of the tree.
   *
   * @param root The root of the tree.
   */
  void Reset(node<node_data> *root);
  /**
   * @brief Get the next node in the traversal.
   *
   * @param history A list to store the history of the traversal.
   * @return A pointer to the next node.
   */
  node<node_data> *Next(_SimpleList *history = NULL);
  /**
   * @brief Get the current node in the traversal.
   *
   * @return A pointer to the current node.
   */
  node<node_data> *Current(void) const;
  /**
   * @brief Get the level of the current node in the tree.
   *
   * @return The level of the current node.
   */
  long Level(void) const;
};

#include "../classes.cp"
#include "../ptr_array.cp"

#endif
