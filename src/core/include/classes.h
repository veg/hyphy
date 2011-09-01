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


//#pragma once
#ifndef __CLASSES__
#define __CLASSES__
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SWAP(a,b,c) ((c)=(a),(a)=(b),(b)=(c))
#define DEBUG 0

#include "hy_lists.h"

/*---------------------------------------------------------------------------*/

template <class array_data> class ptr_array
{

public:
	array_data 	*data;                                            //points to an array of somethings
	int 			length;                                            //length of the array

	ptr_array					(void) {
		data = NULL;
		length = 0;
	}
	//default constructor
	~ptr_array				() {
		if (data) {
			delete [] data;
		}
	}
	//default destructor
	void 			add					(array_data);
	//adds a node to the end of the array
	void	 		prepend				(array_data);
	//adds a node to the beginning of the array
	void 			delete_entry		(int);
	int 			get_length			(void) {
		return length;
	}
	//returns the length of the array  array_data get_elem(int index);
};

/*---------------------------------------------------------------------------*/

template <class node_data> class node
{
public:
	node_data 		in_object;
	ptr_array<node*> 	nodes;
	node 				*parent;
	void 				set_parent		(node &thenode) {
		parent = &thenode;
	}
public:
	node			(void) {
		parent = NULL;
	}
	~node()
	{}
	void 				init			(node_data value) {
		in_object = value;
	}
	node_data 		get_data		(void) {
		return in_object;
	}
	//TEST BELOW
	void 				detach_child	(int k) {
		nodes.delete_entry(k);
	}
#ifndef __UNIX__
	int  				tree_depth 		(void);
#endif
	void 				detach_parent	(void) {
		parent = NULL;
	}
	int  				get_num_nodes	(void) {
		return nodes.get_length();
	}
	void 				add_node		(node& thenode) {
		thenode.set_parent(*this);
		nodes.add(&thenode);
	}

	void				replace_node	(node* existingNode, node* newNode);

	void 				prepend_node	(node &thenode) {
		thenode.set_parent(*this);
		nodes.prepend(&thenode);
	}
	void 				kill_node		(int in) {
		nodes.delete_entry(in);
	}
	node* 			get_node		(int in) {
		return nodes.data[in-1];
	}
	node* 			get_parent		(void) {
		return parent;
	}
	int 				get_child_num	(void);

	node<node_data>* 	go_up			(void);
	node<node_data>*  go_next			(void);
	node<node_data>*  go_previous		(void);
	node<node_data>*  go_down			(int index);
	int 				down			(int index);

	int 				next			(void);
	int 				up				(void);
	void 				delete_tree 	(bool = false);
	node<node_data>*  duplicate_tree  (void);
	// traverse from this node down duplicating the tree and return
	// pointer to the root of the duplicate;
	bool			    compare_subtree (node<node_data>*);

};

template <class node_data>  node<node_data>* StepWiseTraverser 				(node<node_data>* root);
template <class node_data>  node<node_data>* StepWiseTraverserLevel 		(node_data& level, node<node_data>* root);
//template <class array_data> node<node_data>* DepthWiseTraverserLevel 		(node_data& level, node<node_data>* root);
template <class node_data> node<node_data>* DepthWiseStepTraverser 		(node<node_data>* root);
template <class node_data> node<node_data>* DepthWiseStepTraverserRight 	(node<node_data>* root);
template <class node_data> long	    		NodePathTraverser				(_SimpleList& history, node<node_data>* root);
template <class node_data> node<node_data>* DepthWiseStepTraverserWCount  	(long& costCount, node<node_data>* root);
template <class node_data> node<node_data>* DepthWiseStepTraverserLevel    (long& level, node<node_data>* root);

typedef	 union {
	long 		 leafIndex;
	_SimpleList	*leafList;
}
descendantInfo;

node <descendantInfo>*
GatherTreeInfo	 		 		(node<long>* oldRoot, long& leafCounter, _SimpleList& reindex);

void		PurgeTreeInfo					(node<descendantInfo>* root);

#include "../ptr_array.cp"
#include "../classes.cp"

#endif
