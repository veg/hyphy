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

#include "errorfns.h"

template <class node_data>
node<node_data>* StepWiseTraverserLevel (node_data& level, node<node_data>* root)
{
	static node<node_data>* laststep;
	node   <node_data>  * curstep, 
                        * crashdummy;
                        
	bool   goingup = false;
		
	if (root) {
		laststep = root;
		level = 0;
		return root;
	}
	
	curstep = laststep;
	while (curstep) {
		if (!goingup) {
			crashdummy = curstep->go_down(1);
			if (crashdummy) {
				level++;
				curstep=crashdummy;
				break;
			}
		}
		crashdummy = curstep->go_next();
		if (crashdummy) {
			curstep=crashdummy;
			break;
		}
		goingup=true;
		curstep=curstep->go_up();
		level--;
	}
	laststep = curstep;
	return curstep;	
}


//-------------------------------------------------------------

template <class node_data>
node<node_data>* StepWiseTraverser (node<node_data>* root)
{
	long   ignored_level = 0;
	return StepWiseTraverserLevel (ignored_level, root);
}

//-------------------------------------------------------------

template <class node_data>
node<node_data>* DepthWiseStepTraverser  (node<node_data>* root)
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
	curstep=curstep->go_up();
	laststep = curstep;
	return curstep;
}

template <class node_data>
node<node_data>* DepthWiseStepTraverserRight  (node<node_data>* root)
{
	static node<node_data>* laststep;
	node<node_data>* curstep, *crashdummy;
		
	if (root)
	{
		laststep = root;
		while ((crashdummy = laststep->go_down(laststep->get_num_nodes())))
			laststep = crashdummy;
		return laststep;
	}
	
	curstep = laststep;
	crashdummy = curstep->go_previous();
	if (crashdummy)
	{
		curstep=crashdummy;
		while ((crashdummy = curstep->go_down(curstep->get_num_nodes())))
			curstep = crashdummy;
		return laststep = curstep;
	}
	curstep=curstep->go_up();
	laststep = curstep;
	return curstep;
}


template <class node_data>
node<node_data>* DepthWiseStepTraverserWCount  (long& costCount, node<node_data>* root)
{
	static node<node_data>* laststep;
	node<node_data>* curstep, *crashdummy;
		
	if (root)
	{
		laststep = root;
		costCount = root->get_num_nodes();
		while ((crashdummy = laststep->go_down(1)))
		{
			laststep = crashdummy;
			costCount += laststep->get_num_nodes();
		}
		return laststep;
	}
	
	curstep = laststep;
	crashdummy = curstep->go_next();
	if (crashdummy)
	{
		curstep=crashdummy;
		while ((crashdummy = curstep->go_down(1)))
		{
			curstep = crashdummy;
			costCount += curstep->get_num_nodes();
		}
		return laststep = curstep;
	}
	costCount -= curstep->get_num_nodes();
	curstep=curstep->go_up();
	laststep = curstep;
	return curstep;
}
		
		
template <class node_data>
node<node_data>* DepthWiseStepTraverserLevel  (long& level, node<node_data>* root)
{
	static node<node_data>* laststep,* locRoot;
	node<node_data>* curstep, *crashdummy;
		
	if (root)
	{
		laststep = root;
		level = 0;
		while ((crashdummy = laststep->go_down(1))) 
		{
			laststep = crashdummy;
			level++;
		}
		locRoot = root;
		return laststep;
	}
	
	if (laststep==locRoot) return nil;
	
	curstep = laststep;
	crashdummy = curstep->go_next();
	if (crashdummy)
	{
		curstep=crashdummy;
		while ((crashdummy = curstep->go_down(1))) 
		{
			curstep = crashdummy;
			level++;
		}
		return laststep = curstep;
	}
	curstep=curstep->go_up();
	level--;
	laststep = curstep;
	return curstep;
}
		
//-------------------------------------------------------------
/*node <descendantInfo>* 	GatherTreeInfo	 (node<long>* oldRoot, long& leafCounter, _SimpleList& reindex)
{ 		
 	node<descendantInfo>* result = new node<descendantInfo>;
 	
 	checkPointer		 (result);
 	
 	long	nc = oldRoot->get_num_nodes(),
 			i;
 	
 	for (i=1; i<=nc ; i++)
		result->add_node(*GatherTreeInfo(oldRoot->go_down(i),leafCounter,reindex));
 	
 	if (nc)
 	{
 		result->in_object.leafList = new _SimpleList;
 		checkPointer		 (result->in_object.leafList);
 		

 		for (i=1; i<=nc ; i++)
 		{
	 		node<descendantInfo>* thisChild = result->go_down (i);
	 		if (thisChild->get_num_nodes())
	 			(*result->in_object.leafList)<<*(thisChild->in_object.leafList);
	 		else
	 			(*result->in_object.leafList) << thisChild->in_object.leafIndex;
		}
		result->in_object.leafList->Sort();
 	}
 	else
 		result->in_object.leafIndex = reindex.lData[leafCounter++];
 	
 	return result;
}

//-------------------------------------------------------------
bool 	SimpleMatch	 (node<descendantInfo>* root, _SimpleList& toMatch)
{ 		
 	long	nc = root->get_num_nodes(),
 			i;
 	
 	for (i=1; i<=nc ; i++)
 	{
		node<descendantInfo>* aChild = root->go_down (i);
		if (aChild->get_num_nodes())
		{
			_SimpleList * nL = aChild->in_object.leafList;
			if (nL->BinaryFind(toMatch.lData[0])>=0)
			{
				if (nL->lLength < toMatch.lLength)
					return false;
				if (nL->lLength == toMatch.lLength)
					return toMatch.Equal (*nL);
					
				return SimpleMatch (aChild,toMatch);
			}
		}	
	}
 	
  	return false;
}

//-------------------------------------------------------------
void 	PurgeTreeInfo	 (node <descendantInfo>* root)
{ 		
 	long	nc = root->get_num_nodes(),
 			i;
 	
 	for (i=1; i<=nc ; i++)
		PurgeTreeInfo(root->go_down(i));
 	
 	if (nc)
  		DeleteObject (root->in_object.leafList);
  		
  	delete root;
}
*/




















//-------------------------------------------------------------
template <class node_data> void node<node_data>::delete_tree(bool delSelf){
 	
 	 long 	nc = get_num_nodes();
 	 for (int i=1; i<=nc; i++)
 	 {
			go_down(i)->delete_tree(); 	 
			delete (go_down(i));
 	 }
	 if (delSelf)
		 delete (this);
}

//-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::duplicate_tree(){
 		
 	node<node_data>* result = new node<node_data>;
 	for (int i=1; i<=get_num_nodes(); i++) {
		result->add_node(*(go_down(i)->duplicate_tree()));
 	}
 	result->in_object = in_object;
 	return result;
}

//-------------------------------------------------------------
template <class node_data> bool node<node_data>::compare_subtree(node<node_data>* compareTo)
{ 		
	int nNodes = get_num_nodes();
	if (nNodes==compareTo->get_num_nodes())
	{
		for (int i=1; i<=nNodes; i++)
			if (!go_down(i)->compare_subtree (compareTo->go_down(i)))
				return false;
		return true;
	}
	return false; 	
}


//-------------------------------------------------------------

template <class node_data> long NodePathTraverser (_SimpleList& history, node<node_data>* root)
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
}



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
	static int  going_up;
	static node<node_data>* laststep;
	node<node_data>* curstep, *crashdummy;
		
	if (root)
	{
		laststep = root;
		while ((crashdummy = laststep->go_down(1))) laststep = crashdummy;
		going_up = false;
		return laststep;
	}
	
	curstep = laststep;
	crashdummy = curstep->go_next();
	if (crashdummy)
	{
		curstep=crashdummy;
		while ((crashdummy = curstep->go_down(1))) curstep = crashdummy;
		going_up = false;
		return laststep = curstep;
	}
	curstep=curstep->get_parent();
	going_up = true;
	laststep = curstep;
	return curstep;
}
//-----------------------------------Set Number 1----------------


template <class node_data> void node<node_data>::replace_node(node<node_data>* existing, node<node_data>* newNode){
	for (long j = 0; nodes.length; j++)
	{
		if (nodes.data[j] == existing)
		{
			nodes.data[j] = newNode;
			break;
		}
	}	
}

//-----------------------------------Set Number 1----------------
template <class node_data> int node<node_data>::get_child_num()
{
	 int num_siblings;
	 if (parent != NULL)
	   {
		  num_siblings = (*parent).get_num_nodes();
	      for (int i=1; i<num_siblings+1; i++)
	        {
	          if ((*parent).get_node(i) == this) return (i);
	        }
	   }
 	return -1;
}
//----------------------------------end no 1--------------------
//---------Public set no 2-------------------------------------

//--Bool (T/F) responses to potential tree moves
 template <class node_data> int node<node_data>::down(int index)
 {
   if ((index > 0) && (index <= get_num_nodes())){
	  return 1;
   }
   else return 0;
 } //Truth of decent
 
//-------------------------------------------------------------

template <class node_data> int node<node_data>::next(){
   
   if (get_parent() == NULL) return 0;
   if (get_child_num() < (get_parent())->get_num_nodes()) return 1;
   return 0;
 }

//-------------------------------------------------------------
 template <class node_data> int node<node_data>::up(){
   
   if (get_child_num() > 0) return 1;
   else return 0;
 }

//--------MOVERS MOVERS through the tree
template <class node_data> node<node_data>* node<node_data>::go_up(){
  return get_parent();
}

//-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_next(){
  int marker1,marker2;
  marker1 = get_child_num();
  if (get_parent() == NULL) {
    return NULL;
  }
  marker2 = (get_parent())->get_num_nodes();
  if (marker1 < marker2) return ((get_parent())->get_node(marker1+1));
  return NULL;
}

//-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_previous(){
  int marker1;
  if (get_parent() == NULL) {
    return NULL;
  }
  marker1 = get_child_num();
  if (marker1 > 1) return ((get_parent())->get_node(marker1-1));
  return NULL;
}


//-------------------------------------------------------------
template <class node_data> node<node_data>* node<node_data>::go_down(int index){
   if ((index > 0) && (index <= get_num_nodes()))
  	 return (get_node(index));
   
   return NULL; //false=can go no further
 }
//-------------end public set no 2-----------------------------
