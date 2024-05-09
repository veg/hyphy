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

#include <math.h>
#include <ctype.h>
#include <float.h>


#include "global_things.h"
#include "global_object_lists.h"
#include "tree.h"
#include "tree_iterator.h"
#include "hbl_env.h"
#include "category.h"
#include "likefunc.h"

const _String kTreeErrorMessageEmptyTree ("Cannot construct empty trees");


using namespace hy_global;
using namespace hy_env;
using namespace hyphy_global_objects;

#define     TREE_V_SHIFT            8.0
#define     TREE_H_SHIFT            10.0

#ifdef  _SLKP_USE_AVX_INTRINSICS

inline __m256d _hy_matrix_handle_axv_mfma (__m256d c, __m256d a, __m256d b) {
    #if defined _SLKP_USE_FMA3_INTRINSICS
        return _mm256_fmadd_pd (a,b,c);
    #else
        return _mm256_add_pd (c, _mm256_mul_pd (a,b));
    #endif
}

#endif


hyFloat  _TheTree::_timesCharWidths[256]= { // Hardcoded relative widths of all 255 characters in the Times font, for the use of PSTreeString
    0,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0,0.25098,0.721569,0.721569,0.721569,0,0.721569,0.721569,
    0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0,0.721569,0.721569,
    0.25098,0.333333,0.407843,0.501961,0.501961,0.831373,0.776471,0.180392,0.333333,0.333333,0.501961,0.564706,0.25098,0.333333,0.25098,0.278431,
    0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.278431,0.278431,0.564706,0.564706,0.564706,0.443137,
    0.921569,0.721569,0.666667,0.666667,0.721569,0.611765,0.556863,0.721569,0.721569,0.333333,0.388235,0.721569,0.611765,0.890196,0.721569,0.721569,
    0.556863,0.721569,0.666667,0.556863,0.611765,0.721569,0.721569,0.945098,0.721569,0.721569,0.611765,0.333333,0.278431,0.333333,0.470588,0.501961,
    0.333333,0.443137,0.501961,0.443137,0.501961,0.443137,0.333333,0.501961,0.501961,0.278431,0.278431,0.501961,0.278431,0.776471,0.501961,0.501961,
    0.501961,0.501961,0.333333,0.388235,0.278431,0.501961,0.501961,0.721569,0.501961,0.501961,0.443137,0.478431,0.2,0.478431,0.541176,0.721569,
    0.721569,0.721569,0.666667,0.611765,0.721569,0.721569,0.721569,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,
    0.443137,0.443137,0.278431,0.278431,0.278431,0.278431,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,
    0.501961,0.4,0.501961,0.501961,0.501961,0.34902,0.454902,0.501961,0.760784,0.760784,0.980392,0.333333,0.333333,0.54902,0.890196,0.721569,
    0.713725,0.54902,0.54902,0.54902,0.501961,0.576471,0.494118,0.713725,0.823529,0.54902,0.27451,0.27451,0.309804,0.768627,0.666667,0.501961,
    0.443137,0.333333,0.564706,0.54902,0.501961,0.54902,0.611765,0.501961,0.501961,1,0.25098,0.721569,0.721569,0.721569,0.890196,0.721569,
    0.501961,1,0.443137,0.443137,0.333333,0.333333,0.54902,0.494118,0.501961,0.721569,0.168627,0.745098,0.333333,0.333333,0.556863,0.556863,
    0.501961,0.25098,0.333333,0.443137,1,0.721569,0.611765,0.721569,0.611765,0.611765,0.333333,0.333333,0.333333,0.333333,0.721569,0.721569,
    0.788235,0.721569,0.721569,0.721569,0.721569,0.278431,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333
},
_TheTree::_maxTimesCharWidth = 0.980392;

_String const _TheTree::kTreeOutputLabel       ( "TREE_OUTPUT_BRANCH_LABEL"),
              _TheTree::kTreeOutputTLabel      ( "TREE_OUTPUT_BRANCH_TLABEL"),
              _TheTree::kTreeOutputFSPlaceH    ( "__FONT_SIZE__");


#define     DEGREES_PER_RADIAN          57.29577951308232286465


//#define _UBER_VERBOSE_MX_UPDATE_DUMP
//#define _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL 1

extern  long likeFuncEvalCallCount,
              matrix_exp_count;


hyFloat             _lfScalerPower            = 64,
                    _lfScalerUpwards          = pow(2.,_lfScalerPower),
                    _lfScalingFactorThreshold = 1./_lfScalerUpwards,
                    _logLFScaler              = _lfScalerPower *log(2.),
                    _lfMaxScaler              = sqrt(DBL_MAX * 1.e-10),
                    _lfMinScaler              = 1./_lfMaxScaler;



_Vector      _scalerMultipliers,
_scalerDividers;

#ifdef _SLKP_DEBUG_

/*----------------------------------------------------------------------------------------------------------*/
void echoNodeList (_SimpleList& theNodes, _SimpleList& leaves, _SimpleList& iNodes, _SimpleList& flatParents) {
    for (long n = 0; n < theNodes.lLength; n++) {
        bool is_leaf = theNodes(n)<leaves.lLength;
        node<long>* nd = (node<long>*)(is_leaf?leaves(theNodes(n)):iNodes(theNodes(n)-leaves.lLength));
        long children = 0;
        if (!is_leaf) {
            long myIndex = theNodes(n)-leaves.lLength;
            flatParents.Each ([&](long v, unsigned long)->void {if (v == myIndex) children++;});
        }
        printf ("%ld %ld %s (%d)\n", n, theNodes(n), LocateVar(nd->in_object)->GetName()->get_str(), children);
    }
}

#endif

/*----------------------------------------------------------------------------------------------------------*/

__attribute__((__always_inline__)) hyFloat _computeReductionScaler (hyFloat currentScaler, hyFloat sum, long& didScale) {
    if (currentScaler > _lfMinScaler) {
       didScale   = -1;
       sum       *= _lfScalingFactorThreshold;
                                   
       hyFloat tryScale2 = currentScaler * _lfScalingFactorThreshold,
               scaler = _lfScalingFactorThreshold;
       
       while (sum > _lfScalerUpwards && tryScale2 > _lfMinScaler) {
           sum     *= _lfScalingFactorThreshold;
           tryScale2 *= _lfScalingFactorThreshold;
           scaler *= _lfScalingFactorThreshold;
           didScale --;
       }

       return scaler;
    }
    return NAN;
}

__attribute__((__always_inline__)) hyFloat _computeBoostScaler (hyFloat currentScaler, hyFloat sum, long& didScale) {
    if (currentScaler < _lfMaxScaler) {
       didScale                                            = 1;
       sum                                                *= _lfScalerUpwards;
       
       hyFloat tryScale2 = currentScaler * _lfScalerUpwards,
               scaler = _lfScalerUpwards;
        
       while (sum < _lfScalingFactorThreshold && tryScale2 < _lfMaxScaler) {
           sum     *= _lfScalerUpwards;
           tryScale2 *= _lfScalerUpwards;
           scaler *= _lfScalerUpwards;
           didScale ++;
       }
       return scaler;
    } /*else {
#pragma omp critical
        printf ("Overflow in _computeBoostScaler %15.12lg\n", currentScaler);
    }*/
    return NAN;
}
/*----------------------------------------------------------------------------------------------------------*/

hyFloat  acquireScalerMultiplier (long s) {
    if (s>0) {
        if (s >= _scalerMultipliers.get_used())
            for (long k = _scalerMultipliers.get_used(); k <= s; k++) {
                _scalerMultipliers.Store (exp (-_logLFScaler * k));
            }
        return _scalerMultipliers.theData[s];
    }
    s = -s;
    if (s >= _scalerDividers.get_used())
        for (long k = _scalerDividers.get_used(); k <= s; k++) {
            _scalerDividers.Store (exp (_logLFScaler * k));
        }
    return _scalerDividers.theData[s];
}


//_______________________________________________________________________________________________

_TheTree::_TheTree () {
    categoryCount           = 1L;
    aCache                  = nil;
}       // default constructor - doesn't do much



//_______________________________________________________________________________________________

_TheTree::~_TheTree (void) {
      DeleteObject (aCache);

}

//_______________________________________________________________________________________________

void    _TheTree::PurgeTree (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);

    while (_CalcNode * iterator = ti.Next()) {
        DeleteVariable (*iterator->GetName());
    }

    theRoot->delete_tree (true);
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String const & name, _String const & parms, bool make_a_copy):_TreeTopology (&name) {
  PreTreeConstructor   (make_a_copy);
  _TreeTopologyParseSettings settings = CollectParseSettings();
  settings.AllocateCache();
  
  if (_AssociativeList * meta = MainTreeConstructor  (parms, settings)) {
    PostTreeConstructor  (make_a_copy, meta);
  }
}


//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String const & name, _TreeTopology* top):_TreeTopology (&name) {
  PreTreeConstructor   (false);
  if (top->theRoot) {
    isDefiningATree         = kTreeIsBeingParsed;
    theRoot                 = top->theRoot->duplicate_tree ();
    
    _TreeTopologyParseSettings parse_settings = _TreeTopology::CollectParseSettings();
    parse_settings.AllocateCache();
    
    ConditionalTraverser (
        [&] (node<long>* iterator, node_iterator<long> const& ni) -> bool {
          hyFloat   stored_branch_length = top->GetBranchLength (iterator);
          
          _String              nodeVS   ,
                               nodeName (top->GetNodeName    (iterator, false)),
                               *nodeSpec = map_node_to_calcnode(iterator)->GetBranchSpec();

          if (!nodeName.IsValidIdentifier(fIDAllowFirstNumeric)) {
            nodeName  = nodeName.ConvertToAnIdent (fIDAllowFirstNumeric);
          }

          if (stored_branch_length != 0.0) {
            nodeVS = stored_branch_length;
          }
          
          FinalizeNode (iterator, 0, nodeName, *nodeSpec, nodeVS, parse_settings);
          DeleteObject (nodeSpec);
          return false;
        },
        true // SLKP 20180309 : TODO check to see if it is necessary to traverse the root
      );

    isDefiningATree         = kTreeNotBeingDefined;
    PostTreeConstructor      (false, nil);
  } else {
    HandleApplicationError ("Can't create an empty tree");
    return;
  }
}

//_______________________________________________________________________________________________

void    _TheTree::PreTreeConstructor (bool) {
    rooted                  = UNROOTED;
    categoryCount           = 1;
    aCache                  = new _AVLListXL (new _SimpleList);
}

//_______________________________________________________________________________________________

void    _TheTree::delete_associated_calcnode (node<long> * n) const {
   DeleteVariable(*map_node_to_calcnode(n)->GetName());
}


//_______________________________________________________________________________________________

void    _TheTree::PostTreeConstructor (bool make_copy, _AssociativeList* meta) {
  // TODO: SLKP 20180311, extensive duplication with the same function from _TreeToology, consider
  // a possible refactor
  
    if (aCache) {
      DeleteObject (aCache->dataList);
      DeleteAndZeroObject(aCache);
    }
 
  auto variable_handler = [&] (void) -> void {
    /** TODO SLKP 20171211, make sure the semantics are unchanged */
    // existing variable is a CalcNode
    variablePtrs.Replace (get_index(), make_copy ? this->makeDynamic() : this, false);
    setParameter (WrapInNamespace (_TreeTopology::kMeta, GetName()), meta ? meta : new _MathObject, nil, false);
    //printf ("makecopy = %d [%ld]\n", make_copy, this->SingleReference());
  };
  
  bool accept_rooted = EnvVariableTrue(accept_rooted_trees);
  try {
      if (theRoot->get_num_nodes() <= 2) { // rooted tree - check
        if (accept_rooted == false) {
          
          //long node_index = theRoot->get_data();
          bool recurse = false;
          
          if (theRoot->get_num_nodes() == 2) {
            for (int i = 1; i<=2; i++) {
              node<long> *node_temp = theRoot->go_down(i);
              if (node_temp->get_num_nodes()) { // an internal node - make it a root
                delete_associated_calcnode(theRoot);
                node_temp->detach_parent();
                node_temp->add_node(*theRoot->go_down(3-i));
                delete theRoot;
                theRoot = node_temp;
                //delete_associated_calcnode (theRoot);
                rooted = i == 1 ? ROOTED_LEFT : ROOTED_RIGHT;
                ReportWarning (_String("Rooted topology. Removing one branch - the ") & (i==1 ? "left" : "right") & " root child has been promoted to be the new root");
                break;
              }
            }
            
            if (rooted==UNROOTED) {
              ReportWarning ("One branch tree supplied - hopefully this IS what you meant to do.");
              node<long> *delete_me = theRoot->go_down(2);
              //delete_associated_calcnode(theRoot);
              delete_me->detach_parent();
              theRoot->detach_child(2);
              delete_associated_calcnode(delete_me);
              delete delete_me;
              rooted = ROOTED_LEFT;
              ReportWarning ("RIGHT child collapsed to the root");
              //delete_associated_calcnode(theRoot);
            }
          } else {
            if (theRoot->get_num_nodes() == 0) {
              //delete this;
              throw (kTreeErrorMessageEmptyTree);
            }
            node<long> *node_temp = theRoot->go_down(1);
            node_temp->detach_parent();
            delete_associated_calcnode(theRoot);
            delete theRoot;
            theRoot = node_temp;
            ReportWarning ("The root has a single child, which is be promoted to the root");
            recurse = true;
          }
          
          if (recurse) {
            PostTreeConstructor (make_copy, meta);
            return;
          }
        }
      }
      if (!theRoot) {
          //delete this;
          throw _String ("Invalid tree/topology string specification.");
      }
  } catch (const _String& e) {
      HandleApplicationError (e);
  }
    
  variable_handler ();
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String const &name, _TheTree* otherTree):_TreeTopology (&name) {
    PreTreeConstructor   (false);
    if (otherTree->theRoot) {
        isDefiningATree         = kTreeIsBeingParsed;
        theRoot                 = otherTree->theRoot->duplicate_tree ();

        node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        while (node<long>*topTraverser = tree_iterator.Next()) {
            _CalcNode   *sourceNode = (_CalcNode*)LocateVar(topTraverser->in_object),
                          copiedNode (sourceNode, this);
            topTraverser->init (copiedNode.theIndex);
        }

        isDefiningATree         = kTreeNotBeingDefined;
        PostTreeConstructor      (false, nil);
    } else {
        HandleApplicationError(kTreeErrorMessageEmptyTree);
        return;
    }
}

//_______________________________________________________________________________________________

BaseRef  _TheTree::makeDynamic (void) const {
  _TheTree* res = new _TheTree;
  res->_CalcNode::Duplicate (this);
  res->rooted = rooted;
  res->categoryCount = 1;
  res->theRoot = CopyTreeStructure (theRoot,true);
  return res;
}


//_______________________________________________________________________________________________

BaseRef  _TheTree::makeDynamicCopy (_String const* replacementName) const {
  _TheTree* res = new _TheTree;
  
  res->rooted = rooted;
  if (theRoot) {
    _String rn = *replacementName&'.';
    res->theRoot = DuplicateTreeStructure (theRoot, &rn, true);
  } else {
    res->theRoot = nil;
  }
  
  res->SetIndex(variableNames.GetXtra(LocateVarByName (*replacementName)));
  res->theName = new _String (*replacementName);
  res->theName->AddAReference();
  return res;
}

//_______________________________________________________________________________________________

node<long>*  _TheTree::DuplicateTreeStructure (node <long>* source_node_object, _String const* new_name,  bool) const {
  
  // TODO SLKP 20180311: this needs to be reviewed; lots of old ugly code
  
  node<long>* locNode = new node<long>;
  
  for (unsigned long i=1UL; i <= source_node_object->get_num_nodes(); i++) {
    locNode->add_node(*DuplicateTreeStructure (source_node_object->go_down(i), new_name , false));
  }

  _String    replace_me = *GetName()&'.',
             *temp;
  
  _CalcNode* sourceNode = (_CalcNode*)map_node_to_calcnode(source_node_object)->makeDynamic();
  _String    new_node_name = sourceNode->GetName()->Replace(replace_me, *new_name, false);
  _Variable  node_var (new_node_name);
  InsertVar (&node_var);
  locNode->init (node_var.get_index());
  
  sourceNode->ForEachLocalVariable(sourceNode->iVariables,  [&] (long var_idx, long ref_idx, long array_index) {
    _Variable new_node_var (LocateVar (var_idx)->GetName()->Replace(replace_me, *new_name, false));
    (*sourceNode->iVariables)[array_index] = new_node_var.get_index();
  });

  sourceNode->ForEachLocalVariable(sourceNode->dVariables,  [&] (long var_idx, long ref_idx, long array_index) {
    _Variable new_node_var (LocateVar (var_idx)->GetName()->Replace(replace_me, *new_name, false));
    (*sourceNode->dVariables)[array_index] = new_node_var.get_index();
    _Variable * v = LocateVar (new_node_var.get_index());
    //_Formula const *existing_constraint = LocateVar (var_idx)->get_constraint();
    // TODO : this is not a robust way to rename variables
    _String * existing_formula = LocateVar (var_idx)->GetFormulaString(kFormulaStringConversionNormal);
    *existing_formula = existing_formula->Replace (replace_me, *new_name,true);
    _Formula new_constraint (*existing_formula);
    v->SetFormula (new_constraint);
    DeleteObject (existing_formula);
  });

  return locNode;
}

//_______________________________________________________________________________________________

bool _MainTreeConstructor_error (const _String& error, const _String& tree_string, unsigned long index) {
  isDefiningATree = kTreeNotBeingDefined;
  HandleApplicationError (   error & ", in the following string context " &
                tree_string.Cut(index>31L?index-32L:0L,index)&
                " <ERROR HERE> "&
                tree_string.Cut(index+1L,tree_string.length()-index>32L?index+32L:-1L)
             );

  return false;
}


//_______________________________________________________________________________________________


_String    _TheTree::FinalizeNode (node<long>* nodie, long number , _String nodeName, _String const& nodeParameters, _String& nodeValue, _TreeTopologyParseSettings const& settings)
{
    
    static const _String kCommentSuffix ("_comment");
  
    bool isAutoGenerated = nodeName.empty() || (settings.ingore_user_inode_names && nodie->get_num_nodes()>0);
  
    if (isAutoGenerated) {
        nodeName = settings.inode_prefix & number;
    } else {
        if (!nodeName.IsValidIdentifier(fIDAllowFirstNumeric)) {
            if (settings.auto_rename_nodes) {
                _String new_name = nodeName.ConvertToAnIdent(fIDAllowFirstNumeric);
                ReportWarning (_String ("Automatically renamed ") & nodeName.Enquote() & " to " & new_name.Enquote() & " in order to create a valid HyPhy identifier");
                nodeName = new_name;
            } else {
                HandleApplicationError(nodeName.Enquote() & " is not a valid HyPhy identifier (auto-remaing not allowed).");
            }
        }
    }

    _String node_parameters = nodeParameters;

    if (nodie == theRoot) {
        node_parameters = kEmptyString;
        nodeValue       = kEmptyString;
    } else {
        if (node_parameters.empty() ) {
            if (lastMatrixDeclared !=-1L ) {
                node_parameters=* hyphy_global_objects::GetObjectNameByType (HY_BL_MODEL, lastMatrixDeclared, false);
            }
        } else {
            if (settings.parser_namespace.nonempty()) {
                node_parameters = settings.parser_namespace & "." & node_parameters;
            }
        }

        if (node_parameters.nonempty()) {
            ReportWarning (_StringBuffer (128UL) << "Model " << '"' << node_parameters << '"' << " assigned to \"" << nodeName << '"');
        } else {
            ReportWarning (_String("No nodel was assigned to ")& nodeName.Enquote());
        }
    }

    hyTreeDefinitionPhase saveIDT = isDefiningATree;
    isDefiningATree = kTreeNodeBeingCreated;
    _CalcNode cNt (nodeName,node_parameters, 4, this, aCache);
    isDefiningATree = saveIDT;
    nodie->init (cNt.get_index());

    _Constant val (ProcessTreeBranchLength(nodeValue));

    if (nodeValue.nonempty() && settings.accept_user_lengths) {
        if (cNt.CountIndependents () == 1UL) { // can assign default values
          bool use_direct_value = true;
            if (settings.auto_convert_lengths) {
                long nodeModelID = cNt.GetModelIndex();
                if (nodeModelID != HY_NO_MODEL && !IsModelOfExplicitForm(nodeModelID)) {
                    _Formula * expressionToSolveFor = nil;
                    long already_converted = settings.parser_cache->FindLong (nodeModelID);
                  
                    if (already_converted < 0) {
                        _Variable   * tV, * tV2;
                        bool         mByF;
                        RetrieveModelComponents (nodeModelID, tV, tV2, mByF);
                        _String * result = ((_Matrix*)tV->GetValue())->BranchLengthExpression((_Matrix*)tV2->GetValue(),mByF);
                        if (result->nonempty()) {
                            expressionToSolveFor = new _Formula (*result);
                            for (unsigned long cc = 0; cc < cNt.categoryVariables.countitems(); cc++) {
                               _CategoryVariable * thisCC = cNt.get_ith_category (cc);
                              thisCC -> SetValue (new _Constant(thisCC->Mean()), false, true, NULL);
                            }
                            settings.parser_cache->Insert ((BaseRef)nodeModelID, (long)expressionToSolveFor, false, false);
                        }
                        DeleteObject (result);
                    } else {
                        expressionToSolveFor = (_Formula*)settings.parser_cache->GetXtra(already_converted);
                    }

                    if (expressionToSolveFor != nil) {
                        _Variable * solveForMe = LocateVar (cNt.iVariables->list_data[1]);
                        hyFloat modelP = expressionToSolveFor->Brent (solveForMe,solveForMe->GetLowerBound(), solveForMe->GetUpperBound(), 1e-6, nil, val.Value());
                        ReportWarning (_String("Branch parameter of ") & nodeName.Enquote() &" set to " & modelP);
                        cNt.GetIthIndependent(0) ->SetValue(new _Constant (modelP), false,true, NULL);
                        use_direct_value = false;
                    }
                }
            }

            if (use_direct_value) {
                cNt.GetIthIndependent(0)->SetValue (&val,true,true,NULL);
                ReportWarning (_String("Branch parameter of ") & nodeName&" set to " & nodeValue);
            }
        } else {
            ReportWarning (_StringBuffer (128) << nodeName << " has " << _String ((long)(cNt.iVariables?cNt.iVariables->lLength/2:0)) << " parameters - branch length not assigned");
        }
    }

    _CalcNode *nodeVar = (_CalcNode*)LocateVar(cNt.get_index());

    if (nodeVar == NULL) return kEmptyString;

    nodeVar->SetValue (&val,true,true,NULL);

    node_parameters.Clear();
    nodeValue.Clear();

    nodeVar->categoryVariables.TrimMemory();
    nodeVar->categoryIndexVars.TrimMemory();
    nodeVar->_VariableContainer::TrimMemory();

    return nodeName;
}

//_______________________________________________________________________________________________

_String const*    _TheTree::GetNodeModel (node<long>* n) const {
    return map_node_to_calcnode(n)->GetModelName ();
}

//_______________________________________________________________________________________________

const _String    _TheTree::GetNodeName      (node<long>* n, bool fullName) const {
    return GetNodeName(map_node_to_calcnode(n), fullName);
}

//_______________________________________________________________________________________________

const _String    _TheTree::GetNodeName      (_CalcNode const* cn, bool fullName) const {
    if (fullName) {
        return *cn->GetName();
    }
    return cn->GetName()->Cut (GetName()->length()+1L,kStringEnd);
}


//_______________________________________________________________________________________________

BaseRef     _TheTree::toStr (unsigned long) {
  return _TreeTopology::toStr();
}


//__________________________________________________________________________________

void _TheTree::CompileListOfModels (_SimpleList& l) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        long    modelID = iterator->GetModelIndex();
        if (modelID != HY_NO_MODEL && l.Find(modelID) == -1 ) {
            l << modelID;
        }
    }
}

//_______________________________________________________________________________________________
void    _TheTree::SetCompMatrices (long catID) const {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
      iterator->SetCompMatrix (catID);
    }
}

//__________________________________________________________________________________

void _TheTree::SetUp (void) {
  
    flatTree.Clear();
    flatNodes.Clear();
    flatLeaves.Clear();
    flatCLeaves.Clear();
    flatParents.Clear();
  
    _SimpleList flatINodeParents;

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);

    while   (_CalcNode* iterator = ti.Next()) {
        if (ti.IsAtLeaf()) {
          flatCLeaves << iterator;
          flatLeaves << (long)ti.GetNode();
          flatParents << (long)ti.GetNode()->get_parent();
        } else {
            flatTree<<iterator;
            flatNodes<< (long)ti.GetNode();
            flatINodeParents << (long)ti.GetNode()->get_parent();
        }
      
    }

    flatParents << flatINodeParents;

    _SimpleList parentlist (flatNodes),
                indexer (flatNodes.lLength,0,1);

    SortLists   (&parentlist,&indexer);
    flatParents.Each ([&] (long value, unsigned long index) -> void {
      if (value) {
        flatParents[index] = indexer.get(parentlist.BinaryFind(value));
      } else {
        flatParents[index] = -1L;
      }
    });
  
    /*for (unsigned long k=0UL; k<flatParents.countitems(); k++) {
          if (flatParents.get(k)) { // not root
              flatParents[k] = indexer.get(parentlist.BinaryFind(flatParents.get(k)));
          } else {
              flatParents[k] = -1L;
          }
    }*/

 }

//__________________________________________________________________________________

bool _TheTree::AllBranchesHaveModels (long matchSize) const {
  // TODO SLKP 20180313: possible deprecation

  _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);

  while (_CalcNode* iterator = ti.Next()) {
    _Matrix * model = iterator->GetModelMatrix();
    if (! (model && (matchSize < 0L || model->GetHDim()==matchSize))) {
      return false;
    }
  }

  return true;
}


//__________________________________________________________________________________

HBLObjectRef _TheTree::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context, HBLObjectRef cache) {

    switch (opCode) {
       case HY_OP_CODE_TYPE: // Type
        return Type(cache);
        break;
    }

    _MathObject * arg0 = _extract_argument (arguments, 0UL, false);

    if (arg0) {
      switch (opCode) {
        case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
          return TEXTreeString(arg0, cache);
       }

      _MathObject * arg1 = _extract_argument (arguments, 1UL, false);

      if (arg1) {
        switch (opCode) {
          case HY_OP_CODE_PSTREESTRING: //PlainTreeString
            return PlainTreeString(arg0,arg1, cache);
        }
      }
    }

    switch (opCode) {
      case HY_OP_CODE_TEXTREESTRING:
      case HY_OP_CODE_PSTREESTRING:
        WarnWrongNumberOfArguments (this, opCode,context, arguments);
        return new _MathObject;
    }

    return  _TreeTopology::ExecuteSingleOp (opCode,arguments,context);

}




//__________________________________________________________________________________
_String const            _TheTree::GetBranchLengthString (node<long> * n, bool get_expression) const {
    
    _CalcNode* tree_node = map_node_to_calcnode(n);
    
    if (get_expression) {
        
        bool    mbf;

        _Variable *mm,
                *fv;

        RetrieveModelComponents(tree_node->GetModelIndex(), mm, fv, mbf);

        if (mm && fv && mm->ObjectClass() == MATRIX && fv->ObjectClass() == MATRIX) {
            return _String (((_Matrix*)mm->GetValue())->BranchLengthExpression((_Matrix*)fv->GetValue(),mbf));
        } else {
            return kEmptyString;
        }

    } else {
        return tree_node->ComputeBranchLength();
    }
}

//__________________________________________________________________________________
hyFloat            _TheTree::GetBranchLength (node<long> * n) const {
    return map_node_to_calcnode(n)->ComputeBranchLength();
}


//__________________________________________________________________________________
_String const            _TheTree::GetBranchValue (node<long> * n) const {
    hyFloat t = map_node_to_calcnode(n)->Value();
    if (!CheckEqual(t, -1)) {
        return t;
    } else {
        return kEmptyString;
    }
}



//__________________________________________________________________________________
_String const            _TheTree::GetBranchVarValue (node<long> * n, long idx) const {
    _CalcNode * tree_node = map_node_to_calcnode(n);
    
    long iVarValue = tree_node->iVariables->FindStepping(idx,2,1);
    if (iVarValue>0) {
        // model parameter
        return _String(LocateVar (tree_node->iVariables->get(iVarValue-1))->Value());
    } else {
        // user parameter
        _String query = _String('.') & *LocateVar(idx)->GetName();
        
        long local_idx = tree_node->FindLocalVariable (tree_node->iVariables, [&] (long local, long template_var, unsigned long i) -> bool {
            return LocateVar (local)->GetName()->EndsWith(query);
        });
        
        if (local_idx != kNotFound) {
            return tree_node->GetIthIndependent(local_idx)->GetName();
        }
    }
    return kEmptyString;
}


//__________________________________________________________________________________

void    _TheTree::AlignNodes (node<nodeCoord>* theNode) const {
    long k = theNode->get_num_nodes();
    if (k) {
        theNode->in_object.v = (theNode->go_down(1)->in_object.v+theNode->go_down(k)->in_object.v)/2.0;
        theNode->in_object.h = 0;
        for (; k; k--) {
            hyFloat t = theNode->go_down(k)->in_object.h;
            if (t<theNode->in_object.h) {
                theNode->in_object.h = t;
            }
        }
        theNode->in_object.h -= TREE_H_SHIFT;
    } else {
        theNode->in_object.v = 0;
        theNode->in_object.h = 0;
    }
}
//__________________________________________________________________________________

node<nodeCoord>* _TheTree::AlignedTipsMapping (node<long>* iterator, hyFloat& current_offset, bool first, bool respectRoot) const{
    // 20180615 TODO:  SLKP this needs review and possibly deprecation
    if (first) {
        current_offset = 0.0;
        long descendants = theRoot->get_num_nodes();
        node<nodeCoord>* aRoot = new node<nodeCoord>; // reslove rootedness here
        aRoot->in_object.varRef = -1;
        if (rooted == UNROOTED || !respectRoot) {
            for (long k=1L; k<=descendants; k++) {
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(k), current_offset));
            }
            AlignNodes (aRoot);
            return aRoot;
        } else {
            node<nodeCoord>* aChild = new node<nodeCoord>;
            aChild->in_object.varRef = -1;
            if (rooted == ROOTED_LEFT) {
                aRoot->add_node (*aChild);
                for (long k=1L; k<descendants; k++) {
                   aChild->add_node(*AlignedTipsMapping(theRoot->go_down(k), current_offset));
                }
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(descendants), current_offset));
            } else {
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(1), current_offset));
                for (long k=2L; k<=descendants; k++) {
                    aChild->add_node(*AlignedTipsMapping(theRoot->go_down(k), current_offset));
                }
                aRoot->add_node (*aChild);
            }
            AlignNodes (aChild);
            AlignNodes (aRoot);
            return      aRoot;
        }
    } else {
        node<nodeCoord>*aNode = new node<nodeCoord>;
        long descendants = iterator->get_num_nodes();
        for (long k=1L; k<=descendants; k++) {
            aNode->add_node(*AlignedTipsMapping(iterator->go_down(k), current_offset));
        }
        if (!descendants) { // terminal node
            aNode->in_object.v = current_offset;
            aNode->in_object.h = 0;
            current_offset+=TREE_V_SHIFT;
        } else {
            AlignNodes (aNode);
        }
        aNode->in_object.varRef = iterator->in_object;
        return aNode;
    }
}


//__________________________________________________________________________________

void _TheTree::ScaledBranchReMapping (node<nodeCoord>* theNode, hyFloat tw) const
{
    theNode->in_object.h -= tw;
    for (long k=1; k<=theNode->get_num_nodes(); k++) {
        ScaledBranchReMapping (theNode->go_down(k), tw);
    }
}

//__________________________________________________________________________________

hyFloat       _TheTree::DetermineBranchLengthGivenScalingParameter (long varRef, _String& matchString, hyTopologyBranchLengthMode mapMode) const {
    if (mapMode == kTopologyBranchLengthNone) { // was 3
        return 1.;
    }

    _CalcNode * tree_node = (_CalcNode*)LocateVar(varRef);

    hyFloat branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    if (mapMode==kTopologyBranchLengthExpectedSubs) { // was 1
        return tree_node->ComputeBranchLength();
    } else if (mapMode==kTopologyBranchLengthUserLengths) { // was 2
        branchLength = tree_node->Value();
        if (branchLength<=0.0) {
            branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;
        }
    } else {
        auto match_pattern = [&] (long local_idx, long template_index, unsigned long index) -> void {
            _Variable * local_parameter = LocateVar(local_idx);
            if (local_parameter->GetName()->EndsWith(matchString)) {
                branchLength = local_parameter->Compute()->Value();
                if (branchLength <= 0.) {
                    branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;
                }
                throw (0);
            }
        };
        
        try {
            tree_node->ForEachLocalVariable(tree_node->iVariables, match_pattern);
            tree_node->ForEachLocalVariable(tree_node->dVariables, match_pattern);
        } catch (int) {
            return branchLength;
        }
    }
    return branchLength;
}


//__________________________________________________________________________________

node<nodeCoord>* _TheTree::ScaledBranchMapping (node<nodeCoord>* theParent, _String* scalingParameter, long locDepth, long& depth, hyTopologyBranchLengthMode mapMode) const{
    // 20180615 TODO:  SLKP this needs further review and possibly deprecation

    // run a pass of aligned tip mapping then perform one more pass from the root to the children
    // pre-order to remap the length of branches.

    static  hyFloat treeWidth;
    bool     wasRoot = !theParent;

    if (!theParent) { // init
        hyFloat vertical_offset;
        theParent = AlignedTipsMapping (theRoot, vertical_offset, true, true);
        theParent->in_object.h = 0.0;
        treeWidth = 0.;
    }

    node<nodeCoord>* currentN;
    long descendants = theParent->get_num_nodes(),
         k           = 1,
         j,
         b           = -1;

    hyFloat  branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    for  (; k<=descendants; k++) {
        currentN = theParent->go_down(k);
        j        = currentN->in_object.varRef;

        if (j>=0) {

            branchLength  = currentN->in_object.bL = DetermineBranchLengthGivenScalingParameter(j,*scalingParameter,mapMode);
            branchLength += theParent->in_object.h;

            StoreIfGreater(treeWidth, branchLength);
            
            theParent->go_down (k)->in_object.h = branchLength;
            ScaledBranchMapping (theParent->go_down(k), scalingParameter, locDepth+1, depth, mapMode);

        } else {
            theParent->go_down (k)->in_object.h = 0;
            b = k;
        }

    }

    if (k==descendants+1) {
        if (locDepth>=depth) {
            depth = locDepth+1;
        }
    }

    if (wasRoot) {
        if (b>0 && descendants==2) {
            j = (b==1)? 2 : 1;

            ScaledBranchReMapping (theParent->go_down(j),0);
            theParent->go_down(b)->in_object.h = 0;
            ScaledBranchMapping (theParent->go_down(b), scalingParameter, locDepth, depth, mapMode);
        }

        ScaledBranchReMapping (theParent, treeWidth);
        return theParent;
    }
    return nil;
}

//__________________________________________________________________________________

node<nodeCoord>* _TheTree::RadialBranchMapping (node<long>* referenceNode, node<nodeCoord>* parentNode, _String* scalingParameter, hyFloat anglePerTip, long& currentTipID, hyFloat& maxRadius, hyTopologyBranchLengthMode mapMode)
{
    // 20180618 TODO:  SLKP this needs review and possibly deprecation

    // label 1 stores current radial distance from the root
    // label 2 stores the angle of the line to this node
    // h and v store the Cartesian coordinates

    node <nodeCoord>* current_node = new node <nodeCoord>;

    hyFloat          branchL     = 0.,
                        referenceL   = 0.;

    if  (parentNode == nil) {
        current_node->in_object.label1 = 0.0;
        current_node->in_object.label2 = 0.0;
    } else {
        referenceL = parentNode->in_object.label1;
        branchL = DetermineBranchLengthGivenScalingParameter(referenceNode->in_object,*scalingParameter,mapMode);
    }

    long                children    = referenceNode->get_num_nodes();

    current_node->in_object.label1 = referenceL + branchL;
    if (children == 0) {
        current_node->in_object.label2 = anglePerTip * currentTipID++;
    } else {
        hyFloat angleSum = 0.;
        for (long n = 1; n <= children; n++) {
            node<nodeCoord>* newChild = RadialBranchMapping (referenceNode->go_down(n), current_node, scalingParameter, anglePerTip, currentTipID, maxRadius, mapMode);
            current_node->add_node(*newChild);
            angleSum += newChild->in_object.label2;
        }
        current_node->in_object.label2 = angleSum / children;
    }

    current_node->in_object.h        = current_node->in_object.label1*cos(current_node->in_object.label2);
    current_node->in_object.v        = current_node->in_object.label1*sin(current_node->in_object.label2);
    maxRadius = MAX(maxRadius, current_node->in_object.label1);
    current_node->in_object.varRef  = referenceNode->in_object;
    current_node->in_object.bL      = branchL;

    return current_node;
}



//__________________________________________________________________________________

void _TheTree::AssignLabelsToBranches (node<nodeCoord>* theParent, _String* scalingParameter, bool below) {
    // 20180618 TODO:  SLKP this needs review and possibly deprecation

    bool    wasRoot = !(theParent->parent);
    long descendants = theParent->get_num_nodes(),b=-1;

    hyFloat  branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    hyTopologyBranchLengthMode mapMode;
    _String     matchString = DetermineBranchLengthMappingMode(scalingParameter, mapMode);

    for  (long k = 1; k<=descendants; k++) {
        node<nodeCoord>* currentN = theParent->go_down(k);
        long j = currentN->in_object.varRef;
        if (j>=0) {
            branchLength = DetermineBranchLengthGivenScalingParameter(j, matchString, mapMode);
            if (below) {
                currentN->in_object.label2 = branchLength;
            } else {
                currentN->in_object.label1 = branchLength;
            }
            AssignLabelsToBranches (theParent->go_down(k), scalingParameter, below);
        } else {
            if (below) {
                currentN->in_object.label2 = 0.;
            } else {
                currentN->in_object.label1 = 0.;
            }
            b = k;
            AssignLabelsToBranches (theParent->go_down(k), scalingParameter, below);
        }

    }

    if (wasRoot) {
        if ( b>0 && descendants==2 ) {
            
            b = b == 1 ? 2 : 1;
            
            if (below) {
                theParent->in_object.label2 = theParent->go_down(b)->in_object.label2/2.;
                theParent->go_down(b)->in_object.label2/=2.;
            } else {
                theParent->in_object.label1 = theParent->go_down(b)->in_object.label1/2.;
                theParent->go_down(b)->in_object.label1/=2.;
            }
        }
    }
}

//__________________________________________________________________________________

hyFloat _TheTree::PSStringWidth (_String const& s) {

    hyFloat nnWidth = 0.;
    
    s.Each([&] (unsigned char cc, unsigned long) -> void {
        nnWidth += _timesCharWidths [cc];
    });
    
    return nnWidth;
}

//__________________________________________________________________________________
HBLObjectRef _TheTree::PlainTreeString (HBLObjectRef p, HBLObjectRef p2, HBLObjectRef cache) {
    // 20180618 TODO:  SLKP this needs review and possibly deprecation
    
    
  static const _String kTreeOutputEmbed   ("TREE_OUTPUT_EMBED"),
    kTreeOutputOptions          ("TREE_OUTPUT_OPTIONS"),
    kTreeOutputBackground       ( "TREE_OUTPUT_BACKGROUND"),
    kTreeOutputRightMargin      ( "TREE_OUTPUT_RIGHT_MARGIN"),
    kTreeOutputXtraMargin       ( "TREE_OUTPUT_XTRA_MARGIN"),
    kTreeOutputSymbols          ( "TREE_OUTPUT_SYMBOLS"),
    kTreeOutputSymbolSize       ( "TREE_OUTPUT_SYMBOL_SIZE"),
    kTreeOutputExtraPS          ( "TREE_OUTPUT_EXTRA_POSTSCRIPT"),
    kTreeOutputPrefixPS         ( "TREE_OUTPUT_PREFIX_POSTSCRIPT"),
    kTreeOutputLayout           ( "TREE_OUTPUT_LAYOUT");

    
  auto computeChordLength = []  (hyFloat l, hyFloat angle, hyFloat* maxCoord = nil) -> hyFloat {
    hyFloat sinV        = sin(angle),
           cosV         = cos(angle);
    
    if (maxCoord) {
      maxCoord[0] = MAX (maxCoord[0], cosV*l);
      maxCoord[1] = MIN (maxCoord[1], cosV*l);
      maxCoord[2] = MAX (maxCoord[2], sinV*l);
      maxCoord[3] = MIN (maxCoord[3], sinV*l);
    }
    
    return l/MAX(fabs(sinV),fabs(cosV));
  };
  
    _TreeTopologyParseSettings parse_settings = _TreeTopology::CollectParseSettings();
    
    _StringBuffer *     res        = new _StringBuffer (512UL);
    long          treeLayout = 0;

    if (p&&(p->ObjectClass()==STRING)) {
        if (p2&&(p2->ObjectClass()==MATRIX)) {
            node<nodeCoord>* newRoot,
                 *currentNd;

            bool     doEmbed = false;
            bool     doSymbol = false;
            _FString *extraPS   = nil,
                     *prefixPS   = nil;
            long     symbolsize = 3;


            _AssociativeList * toptions  = (_AssociativeList*)FetchObjectFromVariableByType (&kTreeOutputOptions,ASSOCIATIVE_LIST);

            if (toptions) {
                HBLObjectRef lc = toptions->GetByKey (kTreeOutputLayout, NUMBER);
                if (lc) {
                    treeLayout = lc->Value();
                }
                lc = toptions->GetByKey (kTreeOutputEmbed, NUMBER);
                if (lc) {
                    doEmbed = lc->Value();
                }
                lc = toptions->GetByKey(kTreeOutputSymbols, NUMBER);
                if ( lc ) {
                    doSymbol = lc->Value();
                }

                lc = toptions->GetByKey(kTreeOutputExtraPS, STRING);
                if ( lc ) {
                    extraPS = (_FString*)lc->Compute();
                }

                lc = toptions->GetByKey(kTreeOutputPrefixPS, STRING);
                if ( lc ) {
                    prefixPS = (_FString*)lc->Compute();
                }

                lc = toptions->GetByKey(kTreeOutputSymbolSize, NUMBER);
                if ( lc ) {
                    symbolsize = lc->Value();
                }
            }

            _String*        theParam = (_String*) p->toStr(),
                            t;

            bool            scaling         = theParam->length(),
                            doLabelWidth  = true;

            long            tipCount        = 0L,
                            fontSize       = -1L;

            _Matrix*        dimMatrix           = ((_Matrix*)p2->Compute());


            hyFloat      hScale              = 1.0,
                            vScale                = 1.0,
                            labelWidth          = 0.,

                            treeHeight           = (*dimMatrix)(0,1),
                            treeWidth         = (*dimMatrix)(0,0),

                            treeRotation      = (dimMatrix->GetVDim()>2)?(*dimMatrix)(0,2):0.0,
                            treeRadius          ,
                            treeArcAngle        ,
                            totalTreeL          = 0.,

                            mappedTreeHeight = 0.0,
                            mappedTreeHeight2 = 1e+100;

            // letter size in points
            if (treeLayout == 1) {
                treeRadius      = treeWidth;
                treeArcAngle    = treeHeight;

                if (treeRadius <= 0.0) {
                    treeRadius = 300.;
                }
                if (treeArcAngle <= 0.0 || treeArcAngle > 360.) {
                    treeArcAngle = 360.0;
                }

                treeWidth = treeHeight = treeRadius * 2.;

                treeRotation = (treeRotation - floor(treeRotation/360.)*360.)/DEGREES_PER_RADIAN;
            } else {
                if (treeHeight <= 0.0) {
                    treeHeight = 792.;
                }
                if (treeWidth  <= 0.0) {
                    treeWidth  = 576.;
                }
            }



            (*res)<< (_String("%!\n%% PS file for the tree '")&*GetName()&"'.\n%% Generated by "& GetVersionString() & " on " & GetTimeStamp ());
            if (treeLayout == 1) {
                (*res) << "% Radial layout\n"
                 << "/righttext  {dup newpath 0 0 moveto false charpath closepath pathbbox pop exch pop exch sub       4 -1 roll exch sub 3 -1 roll newpath moveto show} def\n";
            }
            if (!doEmbed) {
                (*res)<<"<< /PageSize [" << _String(treeWidth+15) /*wayne added 15 to make trees sit inside the page */
                 <<' '
                 <<_String(treeHeight+15)
                 <<"] >> setpagedevice\n";
            }


            long        xtraChars = 0;

            if (toptions) {
                _Matrix* rgbColor = (_Matrix*)(toptions)->GetByKey (kTreeOutputBackground, MATRIX);
                if (rgbColor) {
                    t = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\nnewpath\n";
                    (*res) << t
                    << "0 0 moveto\n"
                    << "0 "
                    << _String(treeHeight)
                    << " lineto\n"
                    << _String(treeWidth)
                    << ' '
                    << _String(treeHeight)
                    << " lineto\n"
                    << _String(treeWidth)
                    << " 0 lineto\n"
                    << "closepath\nfill\n0 0 0 setrgbcolor\n";
                }

                if ( doSymbol ) {
                    /*add some symbol drawing postscript functions */
                    (*res) << "/size "
                            << _String(symbolsize)
                            << " def\n"
                            << "/box { 0 -0.5 size mul rmoveto\n"
                            << "1 size mul 0 size mul rlineto\n"
                            << "0 size mul 1 size mul rlineto\n"
                            << "-1 size mul 0 size mul rlineto\n"
                            << "closepath\n"
                            << "} def\n"
                            << "/triangle { size size 0.5 mul rlineto 0 size -1 mul rlineto closepath } def\n"
                            << "/circle {currentpoint exch 0.5 size mul add exch 0.5 size mul 180 540 arc\n"
                            << "closepath\n"
                            << "} def\n"
                            << "/diamond { 0 -0.5 size mul rmoveto 0.5 size mul 0 rmoveto 45 rotate 0.707107 size mul 0 rlineto 0 size 0.707107 mul rlineto -0.707107 size mul 0 rlineto -45 rotate  closepath} def\n";
                }


                _Constant* fontSizeIn = (_Constant*)(toptions)->GetByKey (kTreeOutputFSPlaceH, NUMBER);
                if (fontSizeIn) {
                    fontSize = fontSizeIn->Value();
                }

                fontSizeIn = (_Constant*)(toptions)->GetByKey (kTreeOutputRightMargin, NUMBER);
                if (fontSizeIn) {
                    treeWidth = MAX(treeWidth/4+1,treeWidth-fontSizeIn->Value());
                    doLabelWidth = false;
                }

                if ((fontSizeIn = (_Constant*)(toptions)->GetByKey (kTreeOutputXtraMargin, NUMBER))) {
                    xtraChars = fontSizeIn->Value();
                }
            }

            hyTopologyBranchLengthMode    mapMode;
            _String scalingStringMatch;

            if (scaling) {
                scalingStringMatch = DetermineBranchLengthMappingMode (theParam, mapMode);
            }

            if (treeLayout == 1) {
                
                long                tipID = 0;
                EdgeCount (tipCount, tipID);
                tipID               = 0;
                hScale              = 0.;
                newRoot             = RadialBranchMapping (theRoot,nil,&scalingStringMatch,treeArcAngle*pi_const/(180.*tipCount),tipID,hScale,mapMode);
                totalTreeL          = hScale;
            } else {
                if (scaling) {
                    newRoot     = ScaledBranchMapping (nil, &scalingStringMatch, 0, tipCount,mapMode);
                } else {
                    hyFloat offset;
                    newRoot     = AlignedTipsMapping  (theRoot, offset, true);
                }

                hScale    = -treeWidth/newRoot->in_object.h;
            }

            currentNd = NodeTraverser (newRoot);

            while (currentNd) {
                if (currentNd->in_object.v > mappedTreeHeight) {
                    mappedTreeHeight = currentNd->in_object.v;
                }
                if (currentNd->in_object.v < mappedTreeHeight2) {
                    mappedTreeHeight2 = currentNd->in_object.v;
                }

                if (!currentNd->get_num_nodes()) {
                    tipCount++;
                }

                currentNd = NodeTraverser ((node<nodeCoord>*)nil);
            }

            // compute the size of the font, based on the spacing between branches.
            // 36 is max and 2 is min;

            if (fontSize<0) {
                fontSize = treeHeight/(tipCount+1+(scaling>0)*1.5);

                if (fontSize>36) {
                    fontSize = 36;
                } else if (fontSize<2) {
                    fontSize = 2;
                }
            }

            // now recompute the width of the of the tree, including the labels
            // do not allow names to take up more that 1/2 of the width
            // try to keep full length but if the names are too wide,
            // shrink the size of the font

            currentNd = NodeTraverser (newRoot);

            hyFloat  plotBounds[4];

            if (treeLayout == 1) {
                plotBounds [0] = plotBounds[1] = plotBounds [2] = plotBounds[3] = 0.;
            }

            while (currentNd) {
                if (currentNd->in_object.varRef>=0) {
                    _String nodeName (*LocateVar (currentNd->in_object.varRef)->GetName());
                    nodeName.Trim    (nodeName.FindBackwards ('.',0,-1)+1,-1);

                    _AssociativeList * nodeOptions = nil;
                    if (toptions) {
                        nodeOptions = (_AssociativeList *) toptions->GetByKey (nodeName, ASSOCIATIVE_LIST);
                    }


                    HBLObjectRef nodeLabel     = nodeOptions?nodeOptions->GetByKey (_TheTree::kTreeOutputLabel,STRING):nil,
                                nodeTLabel     = nodeOptions?nodeOptions->GetByKey (_TheTree::kTreeOutputTLabel,STRING):nil;

                    hyFloat      nnWidth = 0.0;

                    if (nodeLabel) {
                        nnWidth = 0.0;
                    } else if (nodeTLabel) {
                        nnWidth = 1.+PSStringWidth (((_FString*)nodeTLabel->Compute())->get_str());
                        //printf ("%g\n", nnWidth);
                    } else if (currentNd->get_num_nodes() == 0) {
                        nnWidth = 1.+PSStringWidth (nodeName);
                    }

                    nnWidth += _maxTimesCharWidth * xtraChars;
                    nnWidth *= fontSize;

                    if (treeLayout == 1) {
                        currentNd->in_object.label2 += treeRotation;
                        hyFloat chordLength =  computeChordLength (treeRadius, currentNd->in_object.label2,plotBounds),
                                   overflow    =  MAX(0., treeRadius +
                                                      (nnWidth - chordLength) *
                                                      hScale / MAX(HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS,currentNd->in_object.label1));

                        //nnWidth + currentNd->in_object.label1 * vScale-chordLength);



                        if (overflow > treeRadius*.5) { // shift too big
                            chordLength -= currentNd->in_object.label1/(2.*hScale) * treeRadius;

                            chordLength = chordLength / nnWidth;
                            fontSize = fontSize * chordLength;

                            if (fontSize<2) {
                                fontSize = 2;
                            }

                            nnWidth  = treeRadius*.5;
                        } else {
                            nnWidth = overflow;
                        }
                    } else {
                        if (nnWidth>treeWidth*.5) {
                            fontSize = fontSize / (2.*nnWidth/treeWidth);
                            if (fontSize<2) {
                                fontSize = 2;
                            }
                            nnWidth  = treeWidth*.5;
                        }
                    }
                    if (nnWidth > labelWidth) {
                        labelWidth = nnWidth;
                    }
                }
                currentNd = NodeTraverser ((node<nodeCoord>*)nil);
            }

            if (!doLabelWidth) {
                labelWidth = 0;
            }

            if (scaling) {
                treeHeight -= 3.0*fontSize;
            }

            if (treeLayout == 1) {
                hScale      = (treeRadius-labelWidth)/hScale;
                vScale      = 0.;
            } else {
                hScale      = -(treeWidth-labelWidth)/newRoot->in_object.h;
                vScale      = (treeHeight-2*fontSize)/(mappedTreeHeight - mappedTreeHeight2);
            }
            t = _String ("/Times-Roman findfont\n") & fontSize & " scalefont\nsetfont\n";
            (*res)<<&t;

            nodeCoord dummy;

            long   lw = fontSize/6+1;

            (*res) << _String(lw) << " setlinewidth\n1 setlinecap\n";

             if (prefixPS) {
                (*res) << prefixPS->get_str().Replace (kTreeOutputFSPlaceH, _String(fontSize), true);
            }

            if (treeLayout == 1) {
                //newRoot->in_object.h = -plotBounds[1];
                //newRoot->in_object.v = -plotBounds[3];
                //hScale                 *= 2.*treeRadius / MAX(plotBounds[0]-plotBounds[1],plotBounds[2]-plotBounds[3]);
                newRoot->in_object.h = treeRadius;
                newRoot->in_object.v = treeRadius;
                vScale               = 0.;

                TreePSRecurse (newRoot, (*res), hScale, vScale, ceil(treeWidth), ceil(treeHeight), fontSize/2, labelWidth-fontSize/2, parse_settings, toptions, 1, &treeRadius);
            } else {
                TreePSRecurse (newRoot, (*res), hScale, vScale, ceil(treeWidth), ceil(treeHeight), fontSize/2, labelWidth-fontSize/2, parse_settings, toptions);
            }

            if (scaling) { /* ruler */
                if (fontSize < 8) {
                    fontSize = 8;
                    (*res) << (_String ("/Times-Roman findfont\n") & fontSize & " scalefont\nsetfont\n");
                }

                if (treeWidth < 4*fontSize) { // enforce minimal ruler width
                    treeWidth = 4*fontSize;
                }

                vScale = exp (log(10.)*floor (log10(treeLayout==1?totalTreeL:treeWidth/hScale))); // compute the scaling factor
                if (vScale*hScale < 10.0) {
                    vScale *= 10.;
                }

                _String rulerLabel (vScale, "%5.2g");

                while (vScale*hScale > (treeLayout==1?treeRadius/3:treeWidth-newRoot->in_object.h)) {
                    vScale     *= 0.5;
                    rulerLabel = vScale;
                }

                while (PSStringWidth(rulerLabel)*fontSize>vScale*hScale-3) {
                    vScale      *= 2.0;
                    rulerLabel = vScale;
                }

                long   lm = newRoot->in_object.h,
                       rm = lw+vScale*hScale;

                if (treeLayout == 1) {
                    treeHeight = 2*treeRadius - 2*fontSize;
                    lm = treeWidth - fontSize - rm;
                    rm = treeWidth - fontSize - lw;
                }

                (*res) << "newpath\n"
                      << (_String (lm) & ' ' & _String ((long)treeHeight+2*lw) & " moveto\n")
                      << (_String (rm) & ' ' & (long)(treeHeight+2*lw) & " lineto\nstroke\nnewpath\n") // main horizontal bar
                      << (_String (lm) & ' ' & _String ((long)treeHeight) & " moveto\n")
                      << (_String (lm) & ' ' & (long)(treeHeight+4*lw) & " lineto\nstroke\nnewpath\n")// left notch
                      << ( _String (rm) & ' ' & (long)(treeHeight) & " moveto\n")
                      << (_String (rm) & ' ' & (long)(treeHeight+4*lw) & " lineto\nstroke\nnewpath\n") // right notch
                      << (_String (lm+lw) & ' ' & _String ((long)treeHeight+3*lw) & " moveto\n") // main text
                      << (_String ('(') & _String (rulerLabel) & ") show\n");
            }

            newRoot->delete_tree ();
            delete  newRoot;

            if (extraPS) {
                (*res) << extraPS->get_str().Replace (kTreeOutputFSPlaceH, _String(fontSize), true);
            }

            if (!doEmbed) {
                (*res)<< "showpage";
            }
            DeleteObject (theParam);
        } else {
            ReportWarning ("An invalid 3rd parameter was passed to PSTreeString.");
        }
    } else {
        ReportWarning ("An invalid 2nd parameter was passed to PSTreeString.");
    }
    res->TrimSpace();
    return _returnStringOrUseCache(res, cache);
}




//_______________________________________________________________________________________________

long    _TheTree::CountTreeCategories (void) {
    // 20180618 TODO:  SLKP this is an expensive call (rescans all category variables recursively);
    // could use an optimization pass

    categoryVariables.Clear();
    _AVLList           cVA (&categoryVariables);
    ScanForCVariables (cVA);
    cVA.ReorderList   ();
    categoryCount = 1L;
    
    categoryVariables.Each ([&] (long v, unsigned long) -> void {
        categoryCount *= ((_CategoryVariable*)LocateVar(v))->GetNumberOfIntervals();
    });
    return categoryCount;
}


//__________________________________________________________________________________

void    _TheTree::TreePSRecurse (node<nodeCoord>* iterator, _StringBuffer &res, hyFloat hScale, hyFloat vScale,
                                 long hSize, long vSize, long halfFontSize, long shift, const _TreeTopologyParseSettings& settings, _AssociativeList* outOptions,
                                 char layout, hyFloat * xtra) const {
    
    // TODO: SLKP 20180803, this needs review and possible deprecation
    
    
    static _String const kTreeOutputThickness        ( "TREE_OUTPUT_BRANCH_THICKNESS"),
                         kTreeOutputLinecap          ( "TREE_OUTPUT_BRANCH_LINECAP"),
                         kTreeOutputSplit            ( "TREE_OUTPUT_BRANCH_SPLIT"),
                         kTreeOutputNotchesColor     ( "TREE_OUTPUT_BRANCH_NOTCHES_COLOR"),
                         kTreeOutputNotches          ( "TREE_OUTPUT_BRANCH_NOTCHES"),
                         kTreeOutputColor            ( "TREE_OUTPUT_BRANCH_COLOR"),
                         kTreeOutputDash             ( "TREE_OUTPUT_BRANCH_DASH"),
                         kTreeOutputOLabel           ( "TREE_OUTPUT_OVER_BRANCH"),
                         kTreeOutputNNPlaceH         ( "__NODE_NAME__");

    

    
    // 20180618 TODO:  SLKP this needs review and possibly deprecation

    unsigned long               descendants = iterator->get_num_nodes();
    bool                        is_leaf     = descendants == 0UL;
    //lineW    = halfFontSize/3+1;

    hyFloat         vc,
                       hc,
                       vcl,
                       hcl,
                       hc1,
                       hc2;

    _String            t,
                       varName,
                       colorString ("0 0 0 setrgbcolor\n");

    if (!iterator->is_root()) {
        if (iterator->in_object.varRef >=0) {
          varName = *map_node_to_calcnode(iterator)->GetName();
        }
    } else {
        hc = GetRoot().in_object;
        if (hc >= 0.) {
            varName = (*LocateVar(hc)->GetName());
        }
        if (layout == 1) {
            res <<  (_String (iterator->in_object.h) & ' ' & _String (iterator->in_object.v) & " translate\n");
        }
    }

    varName.Trim(varName.FindBackwards ('.',0L,-1L)+1L,-1L);

    _AssociativeList * nodeOptions = outOptions ? (_AssociativeList *) outOptions->GetByKey (varName, ASSOCIATIVE_LIST) : nil;

    HBLObjectRef nodeLabel  = nodeOptions?nodeOptions->GetByKey (_TheTree::kTreeOutputLabel,STRING):nil,
                 nodeTLabel = nodeOptions?nodeOptions->GetByKey (_TheTree::kTreeOutputTLabel,STRING):nil;

    if (layout == 1) {
        hcl = iterator->in_object.h*hScale;
        vcl = iterator->in_object.v*hScale;
    } else {
        vcl = vSize-iterator->in_object.v*vScale,
        hcl = hSize+iterator->in_object.h*hScale-shift;
    }

    if (is_leaf || nodeLabel)
        // terminal node or default label
    {
        t = kEmptyString;
        hyFloat myAngle = layout==1?iterator->in_object.label2*DEGREES_PER_RADIAN:0.0;
        if (layout == 1) {
            res << (_String(myAngle) & " rotate\n");
            vc = 0;
            hcl = (vScale + iterator->in_object.bL)*hScale;
            hc  = hcl + halfFontSize;
        } else {
            vc = vcl-halfFontSize;
            hc = hcl+halfFontSize;
        }

        if (!nodeTLabel) {
            res<<"newpath\n";

            if (nodeLabel) {
                res << (_String(hcl) & ' ' & _String (vc) & " moveto\n");
                t = ((_FString*)nodeLabel->Compute())->get_str().Replace (kTreeOutputNNPlaceH, varName, true).Replace (kTreeOutputFSPlaceH, _String(halfFontSize*2), true) & '\n';
            }

            if (is_leaf && t.empty()) {
                // generate the default label
                if (layout == 1 && myAngle > 90. && myAngle < 270.) {
                    hyFloat xt = hc-halfFontSize/2,
                               yt = vc-2*halfFontSize/3;
                    t = _String (xt) & _String (" 0 translate 180 rotate 0 ") & _String (yt) & _String ('(') & varName & ") righttext 180 rotate -" & xt & " 0 translate\n";
                } else {
                    res<< (_String(hc-halfFontSize/2) & ' ' & _String (vc-2*halfFontSize/3) & " moveto\n");
                    t = _String ('(') & varName & ") show\n";
                }
            }
            res<<&t;
        }

        if (is_leaf) {
            iterator->in_object.h = hc-halfFontSize;
        }

        if (layout == 1) {
            res << (_String(-iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
        }
    }

    long  minChildHC = 0x0fffffff,
          newV       = 0;

    if (!is_leaf) {
        vc = vSize - iterator->in_object.v*vScale;
        hc = hSize + iterator->in_object.h*hScale-shift;

        nodeCoord childCoord;
        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = iterator->go_down(k);
            TreePSRecurse (child, res, hScale, (layout==1)?vScale+iterator->in_object.bL:vScale, hSize, vSize,halfFontSize,shift,settings, outOptions,layout,xtra);
            if (k==1) {
                hc1 = layout==1?child->in_object.label2:child->in_object.v;
            }
            if (k==descendants) {
                hc2 = layout==1?child->in_object.label2:child->in_object.v;
            }
        }

        char doVLines = 3;

        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = iterator->go_down(k);

            if (child->in_object.varRef>=0) {
                t = map_node_to_calcnode(child)->ContextFreeName();
            } else {
                t = kEmptyString;
            }

            newV += child->in_object.v;

            _AssociativeList * childOptions = nil;

            hyFloat         splitBranch = -1.,
                               lineWP = 0.0;

            _String            childColor,
                               notchColor,
                                childDash,
                                linewidth1,
                                linewidth2,
                                linecap1,
                                linecap2;
            
            _String const *     blabelString = nil;


            _Matrix         *  notches    = nil,
                               *  multiColor = nil;

            if (layout == 1) {
                res << (_String(child->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }

            if (outOptions) {
                childOptions = (_AssociativeList *) outOptions->GetByKey (t, ASSOCIATIVE_LIST);
                if (childOptions) {
                    HBLObjectRef keyVal = childOptions->GetByKey (kTreeOutputThickness,NUMBER);
                    if (keyVal) {
                        lineWP = keyVal->Compute()->Value();
                        linewidth1 = _String("currentlinewidth ") & lineWP & " setlinewidth\n";
                        linewidth2 = "setlinewidth\n";
                    }
                    keyVal = childOptions->GetByKey (kTreeOutputLinecap,NUMBER);
                    if (keyVal) {
                        linecap1 = _String("currentlinecap ") & (long)keyVal->Compute()->Value() & " setlinecap\n";
                        linecap2 = "setlinecap\n";
                    }
                    keyVal = childOptions->GetByKey (kTreeOutputSplit,NUMBER);
                    if (keyVal) {
                        splitBranch = keyVal->Compute()->Value();
                    }

                    keyVal = childOptions->GetByKey (kTreeOutputNotches,MATRIX);
                    if (keyVal) {
                        notches = (_Matrix*)(((_Matrix*)keyVal->Compute())->ComputeNumeric())->makeDynamic();
                    }

                    keyVal = childOptions->GetByKey (kTreeOutputColor,MATRIX);
                    if (keyVal) {
                        _Matrix* rgbColor = (_Matrix*)keyVal->Compute();
                        if (rgbColor->GetHDim() > 1 && rgbColor->GetVDim () == 4 && layout != 1) {
                            multiColor = (_Matrix*)rgbColor->makeDynamic();
                        } else {
                            childColor = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\n";
                        }
                    }

                    keyVal = childOptions->GetByKey (kTreeOutputNotchesColor,MATRIX);
                    if (keyVal) {
                        _Matrix* rgbColor = (_Matrix*)keyVal->Compute();
                        notchColor = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\n";
                    }

                    keyVal = childOptions->GetByKey (kTreeOutputDash,MATRIX);
                    if (keyVal) {
                        _Matrix* dash = (_Matrix*)keyVal->Compute();
                        childDash = _String('[') & _String((*dash)(0,0)) & " " & _String((*dash)(0,1)) & "] " & _String((*dash)(0,2)) & " setdash\n";
                    }
                    keyVal = childOptions->GetByKey (kTreeOutputOLabel,STRING);
                    if (keyVal) {
                        blabelString = &((_FString*)keyVal->Compute())->get_str();
                    }
                }
            }


            if (blabelString) {
                if (layout == 1) {
                    t = _String(iterator->in_object.label1*hScale) & " 0 moveto\n";
                } else {
                    t = _String(hc) & ' ' & _String (child->in_object.v) & " moveto\n";
                }
                res<<&t
                    << (blabelString->Replace (kTreeOutputNNPlaceH, varName, true).Replace (_TheTree::kTreeOutputFSPlaceH, _String(halfFontSize*2), true)) << '\n'
                    <<'\n';
            }

            res << childDash << childColor << linewidth1 << linecap1;

            if (layout == 1) {
                if (child->get_num_nodes()) { // internal node
                    res << "newpath\n" << ((_String("0 0 ") & child->in_object.label1*hScale & ' ' & child->in_object.v & ' ' & (child->in_object.auxD) & " arc \n")) << "stroke\n";
                }
            } else {
                if (childOptions && descendants == 2) {
                    res << "newpath\n"
                        << ((_String(hc) & ' ' & _String (0.5*(hc1+hc2)) & " moveto\n"))
                        << ((_String(hc) & ' ' & _String (k==1?hc1:hc2) & " lineto\n"))
                        << "stroke\n";
                    
                    doVLines -= k;
                }

            }

            res << "newpath\n";
            if (layout == 1) {
                res << (_String(child->in_object.label1*hScale) & " 0 moveto\n");
                t = _String(-child->in_object.bL*hScale) & " 0 rlineto\n";
            } else {

                hyFloat lineWidthInset = 0.0;

                //if (lineWP > 0.0)
                //  lineWidthInset = (lineWP-lineW)*.5;

                if (multiColor) {
                    hyFloat span        = child->in_object.h - hc + 2*lineWidthInset,
                               currentX    = hc - lineWidthInset;

                    res <<  (_String(currentX) & ' ' & _String (child->in_object.v) & " moveto\n");
                    for (long seg = 0; seg < multiColor->GetHDim(); seg++) {
                        res << (_String((*multiColor)(seg,0)) & " " & _String((*multiColor)(seg,1)) & " " & _String((*multiColor)(seg,2)) & " setrgbcolor\n");
                        hyFloat mySpan = span*(*multiColor)(seg,3);
                        res << (_String(mySpan) & " 0 rlineto\n");
                        if (seg < multiColor->GetHDim()-1) {
                            currentX += mySpan;
                            res << "stroke\nnewpath\n"
                                <<  ((_String(currentX) & ' ' & _String (child->in_object.v) & " moveto\n"));
                        }
                    }
                    DeleteObject (multiColor);
                    t = kEmptyString;
                } else {
                    res<< (_String(hc - lineWidthInset) & ' ' & _String (child->in_object.v) & " moveto\n");
                    t = _String(child->in_object.h) & ' ' & _String (child->in_object.v + lineWidthInset) & " lineto\n";
                }
                if (layout == 1) {
                    minChildHC = MIN(child->in_object.label1,minChildHC);
                } else {
                    minChildHC = MIN(child->in_object.h,minChildHC);
                }
            }
            res <<t
                << "stroke\n"
                << linecap2
                << linewidth2;

            if (childDash.nonempty()) {
                res << "[] 0 setdash\n";
            }

            if (splitBranch >= 0.0 && splitBranch <= 1.0) {
                res << "newpath\n";
                hyFloat x,
                           y;

                if (layout == 1) {
                    x = (child->in_object.label1 - child->in_object.bL*splitBranch)*hScale;
                    y = 0.;
                } else {
                    x = hc+(child->in_object.h-hc)*(1.-splitBranch);
                    y = child->in_object.v;
                }

                res<< (_String(x) & ' ' & _String (y) & " " & halfFontSize  & " 0 360 arc\n") << "fill\n";
            }

            if (notches) {
                notches->CheckIfSparseEnough(true);
                res << notchColor;
                for (long l = 0; l < notches->GetSize(); l++) {
                    hyFloat aNotch = (*notches)[l];
                    if (aNotch >= 0. && aNotch <= 1.) {
                        res << "newpath\n";
                        hyFloat x,
                                   y;

                        if (layout == 1) {
                            x = (child->in_object.label1 - child->in_object.bL*aNotch)*hScale;
                            y = 0.;
                        } else {
                            x = hc+(child->in_object.h-hc)*(1.-aNotch);
                            y = child->in_object.v;
                        }

                        res << (_String(x-0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " moveto ")
                         << (_String(x+0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " lineto\n")
                         << (_String(x-0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " moveto ")
                         << (_String(x+0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " lineto\n")
                         << "stroke\n";
                    }
                }
                DeleteObject (notches);
            }

            res << &colorString;
            if (layout == 1) {
                res << (_String(-child->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }
        }

        if (layout == 0 && doVLines) {
            _String linewidth1,
                    linewidth2,
                    linecap1,
                    linecap2;

            if (nodeOptions) {
                HBLObjectRef keyVal = nodeOptions->GetByKey (kTreeOutputThickness,NUMBER);
                if (keyVal) {
                    hyFloat lineWP = keyVal->Compute()->Value();
                    linewidth1 = _String("currentlinewidth ") & lineWP & " setlinewidth\n";
                    linewidth2 = "setlinewidth\n";
                }
                keyVal = nodeOptions->GetByKey (kTreeOutputLinecap,NUMBER);
                if (keyVal) {
                    linecap1 = _String("currentlinecap ") & (long)keyVal->Compute()->Value() & " setlinecap\n";
                    linecap2 = "setlinecap\n";
                }
            }

            res << linecap1
                << linewidth1
                << "newpath\n";

            if (doVLines == 3) {
                t = _String(hc) & ' ' & _String (hc2) & " moveto\n";
            } else {
                t = _String(hc) & ' ' & _String (0.5*(hc1+hc2)) & " moveto\n";
            }
            res<< t;
            if (doVLines == 3) {
                t = _String(hc) & ' ' & _String (hc1) & " lineto\n";
            } else {
                t = _String(hc) & ' ' & _String (doVLines==1?hc1:hc2) & " lineto\n";
            }
            res<< t
             << "stroke\n"
             << &linewidth2
             << &linecap2;
        }

        if (layout == 0) {
            iterator->in_object.h = hc;
            iterator->in_object.v = newV/descendants;
        } else {
            iterator->in_object.auxD = (hc2-iterator->in_object.label2)*DEGREES_PER_RADIAN;
            if (!(iterator->is_root())) {
                iterator->in_object.v    = (hc1-iterator->in_object.label2)*DEGREES_PER_RADIAN;
            }
        }
    } else {
        iterator->in_object.v = vc;
    }

    if (nodeTLabel) {
        t = ((_FString*)nodeTLabel->Compute())->get_str();
        if (t.nonempty()) {
            hyFloat    scF     = 2.*halfFontSize;

            if (layout == 1) {
                res << (_String(iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            } else {
                if (!is_leaf) {
                    vcl     = iterator->in_object.v;
                }
            }


            if (scF != 2.*halfFontSize) {
                res << (_String("/Times-Roman findfont\n")& scF & " scalefont\nsetfont\n");
            }

            res << "newpath\n";
            if (layout == 1) {
                res << (_String(iterator->in_object.label1 * hScale + halfFontSize) & ' ' & _String (-2*halfFontSize/3) & " moveto\n");
            } else {
                if (!is_leaf) {
                    hc = hcl + scF*0.5;
                    vc = vcl - scF*0.33;

                } else {
                    hc -= scF*0.25;
                    vc -= scF*0.33;
                }
                res << (_String(hc) & ' ' & _String (vc) & " moveto\n");
            }

            res << '(' << t << ") show \n";

            if (scF != 2.*halfFontSize) {
                res << (_String("/Times-Roman findfont\n")& halfFontSize*2 & " scalefont\nsetfont\n");
            }

            if (layout == 1) {
                res << (_String(-iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }
        }
    }

    if (colorString.nonempty()) {
        res << "0 0 0 setrgbcolor\n";
    }

    if (iterator->is_root() == false && layout == 1) {
        res <<  (_String (-iterator->in_object.h) & ' ' & _String (-iterator->in_object.v) & " translate\n");
    }
}

//__________________________________________________________________________________

void _TheTree::SetUpMatrices (long categCount) {
  //fprintf (stderr, "[_TheTree::SetUpMatrices] %ld\n", categCount);

    categoryCount = Maximum (categCount,1L);

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);

    while   (_CalcNode* iterator = ti.Next()) {
        if (iterator->IsConstant()) {
            iterator->varFlags |= HY_VC_NO_CHECK;
        }
        iterator->ConvertToSimpleMatrix(categoryCount);

        if (categoryCount==1L) {
            iterator->matrixCache = nil;
        } else {
            iterator->matrixCache = new _Matrix* [categoryCount];
            InitializeArray (iterator->matrixCache, categCount, (_Matrix*)nil);
        }
    }
}


//__________________________________________________________________________________

void _TheTree::CleanUpMatrices (long category_count) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);

    if (categoryCount == 1L) {
        while   (_CalcNode* iterator = ti.Next()) {

            // mod 05/03/2003 - uncomment next 5 lines
            // this breaks after ReplicateConstraint or MolecularClock is called
            // WTF?

            iterator->ConvertFromSimpleMatrix(categoryCount);

            if (iterator->compExp) {
                DeleteObject (iterator->compExp);
                iterator->compExp = nil;
            }

            iterator->varFlags &= HY_VC_CLR_NO_CHECK;
        }
    } else {
        while   (_CalcNode* iterator = ti.Next()) {
            iterator->ConvertFromSimpleMatrix(category_count);

            for (long i=0; i<categoryCount; i++) {
                DeleteAndZeroObject(iterator->matrixCache[i]);
            }

            delete [] iterator->matrixCache;
            iterator->matrixCache = nil;
            iterator->compExp = nil;
            iterator->varFlags &= HY_VC_CLR_NO_CHECK;
        }
        categoryCount = 1;
    }
}

//__________________________________________________________________________________

void _TheTree::RemoveModel (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        iterator->RemoveModel();
    }
    categoryCount = 1;
}



//__________________________________________________________________________________

void _TheTree::ScanContainerForVariables (_AVLList& l,_AVLList& l2, _AVLListX * tagger, long weight, _AVLListX * map_variables_to_nodes) const {
    /**
        map_variables_to_nodes, if supplied, will be used to map variable IDs to "post-order" node index of nodes that they affect, IF they only affect one node; otherwise they will be tagged with '-1'
     */
    
    unsigned long traversal_order = 0UL,
                  leaf_index = 0UL,
                  int_index = 0UL;
    
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
        bool is_leaf = ti.IsAtLeaf();
        
        
        iterator->ScanContainerForVariables(l,l2, tagger, weight +  flatNodes.lLength + flatLeaves.lLength - traversal_order, map_variables_to_nodes, is_leaf ? leaf_index : int_index + flatLeaves.lLength);
        
        if (is_leaf) {
            leaf_index ++;
        } else {
            int_index ++;
        }
        traversal_order ++;
    }
}

//__________________________________________________________________________________

void _TheTree::ScanAndAttachVariables (void) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
        iterator->ScanAndAttachVariables();
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForDVariables (_AVLList& l,_AVLList& l2) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
      iterator->ScanForDVariables(l,l2);
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForGVariables (_AVLList& li, _AVLList& ld, _AVLListX * tagger, long weight) const {
    _SimpleList cL;
    _AVLList    cLL (&cL);
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot );
    
     
    while   (_CalcNode* iterator = ti.Next()) {

        _Formula *explicitFormMExp = iterator->GetExplicitFormModel ();
        _Matrix  *modelM = explicitFormMExp?nil:iterator->GetModelMatrix();

        if ((explicitFormMExp && cLL.Find ((BaseRef)explicitFormMExp) < 0) || (modelM && cLL.Find(modelM) < 0)) {
            _SimpleList temp;
            _AVLList tempA (&temp);
            if (modelM) {
                modelM->ScanForVariables(tempA, true);
                tempA.ReorderList();
                
                for (unsigned long i=0; i<temp.lLength; i++) {
                    long p = temp.list_data[i];
                    _Variable* v = LocateVar (p);
                    if (v&&v->IsGlobal()) {
                        if(v->IsIndependent()) {
                            li.Insert ((BaseRef)p);
                            if (tagger) {
                                tagger->UpdateValue((BaseRef)p, weight, 1);
                            }
                        } else {
                            ld.Insert ((BaseRef)p);
                        }
                    }
                }
            } else {
                explicitFormMExp->ScanFForVariables(tempA, true, false, true, true, tagger, weight, weight >> 2);
            }
            
            cLL.Insert (modelM?(BaseRef)modelM:(BaseRef)explicitFormMExp);
        }
        iterator -> ScanForGVariables(li,ld);
     }
}

//__________________________________________________________________________________

void _TheTree::ScanForCVariables (_AVLList& lcat) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
        for (unsigned long i = 0UL; i < iterator->categoryVariables.lLength; i++) {
            lcat.Insert ((BaseRef)iterator->categoryVariables.Element(i));
        }
    }
}

//__________________________________________________________________________________

bool _TheTree::HasChanged (bool, _AVLListX*) {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
    while   (_CalcNode* iterator = ti.Next()) {
        if (iterator->HasChanged()) {
            return true;
        }
    }
    return false;
}

  //__________________________________________________________________________________

bool _TheTree::HasChanged2 (void) {

  for (unsigned long k = 0; k < categoryVariables.lLength;  k++) {
    if (((_CategoryVariable*)LocateVar(categoryVariables.Element(k)))->HaveParametersChanged()) {
      return true;
    }
  }

  _SimpleList _c;
  _AVLListX cache (&_c);
    
  _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);
  while   (_CalcNode* iterator = ti.Next()) {
    if (iterator->_VariableContainer::HasChanged(false, &cache)) {
      return true;
    }
  }
  return false;
}


//_______________________________________________________________________________________________
void     _TheTree::InitializeTreeFrequencies (_Matrix *mx, bool setDim) {
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the pruning algo
// must be called before any tree pruning computations are started

    unsigned long vecDim = mx->GetHDim()*mx->GetVDim();
    // theModel = mx;

    if  (setDim) {
        SetTreeCodeBase (vecDim);
    } else {
        CopyArray (theProbs, mx->theData, vecDim);
    }
}

//_______________________________________________________________________________________________
void     _TheTree::SetTreeCodeBase (long b) {
  // this will take the  matrix of frequencies and
  // 1) use its dimensions to initialize tree freq holders
  // 2) place global frequencies into the tree holder for later use by the pruning algo
  // must be called before any tree pruning computations are started
  SetCodeBase (b);
  if (cBase>0) {
      _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
      while   (_CalcNode* iterator = ti.Next()) {
        iterator -> SetCodeBase (b);
      }
  }
}

//_______________________________________________________________________________________________

long     _TheTree::IsLinkedToALF (long& pid) const {
    return likeFuncList.FindOnCondition([&] (BaseRefConst value, unsigned long index) -> bool {
        if (likeFuncNamesList.GetItem(index)) {
            pid = ((_LikelihoodFunction const*)value)->DependOnTree (*GetName());
            if (pid >= 0) {
                return true;
            }
        }
        return false;
    });
}



//_______________________________________________________________________________________________

bool     _TheTree::IntPopulateLeaves (_DataSetFilter const* dsf, long site_index) const {
// assign proper values to leaf conditional probability vectors
    bool site_has_all_gaps = true;

    _String* buffer = dsf->MakeSiteBuffer();

    for (long leaf_index = 0; leaf_index<flatLeaves.lLength; leaf_index++) {
        _CalcNode * iterator = (_CalcNode*)flatCLeaves.GetItem(leaf_index);
        dsf->RetrieveState(site_index, leaf_index, *buffer, false);

      site_has_all_gaps &= (dsf->Translate2Frequencies (*buffer, iterator->theProbs, true)<0); // ambig
      site_has_all_gaps &= (!ArrayAny (iterator->theProbs, cBase, [](hyFloat x, unsigned long) {return x == 0.0; })); //completely unresolved
      map_node_to_calcnode (((node <long>*)flatLeaves.GetElement(leaf_index))->parent)->cBase = -1;
    }

    DeleteObject (buffer);
    return site_has_all_gaps;
}


//_______________________________________________________________________________________________

void     _TheTree::RecoverNodeSupportStates (_DataSetFilter const* dsf, long site_index, _Matrix& resultMatrix) {
// TODO SLKP 20180803 this needs to be moved to using standard log likelihood calculations

// assume current values of all parameters
// return 2 sets of vectors for each branch
//   - top-down  conditionals
//   - bottom-up conditionals
//   resultMatrix is assumed to contain
//      uniqueSites X (flatLeaves.lLength+flatTree.lLength)*cBase*2 X categoryCount

    long      globalShifter        = (flatLeaves.countitems()+flatTree.countitems())*cBase,
              catShifer             = dsf->GetPatternCount() * 2 * globalShifter;
    
    IntPopulateLeaves (dsf, site_index);
    
    /* pass 1; populate top-down vectors */
    /* ugly top-bottom algorithm for debuggability and compactness */
    
    // SLKP 20180803 -- this is a very ugly hack to avoid reitroducing CalcNode->nodeIndex
    // this should be taken care of when RecoverNodeSupportStates is re-engineered
    _SimpleList _avl_storage;
    _AVLListX    node_to_index (&_avl_storage);
    long leaf_count = 0L, int_node_count = flatLeaves.countitems();
    
    node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>*topTraverser = tree_iterator.Next()) {
        _CalcNode   *sourceNode = map_node_to_calcnode (topTraverser);
        if (topTraverser->is_leaf()) {
            node_to_index.Insert (sourceNode, leaf_count++, false, false);
        } else {
            node_to_index.Insert (sourceNode, int_node_count++, false, false);
        }
    }
    
    
    for (long catCount   = 0L; catCount < categoryCount; catCount ++) {
        hyFloat * currentStateVector = resultMatrix.theData + 2*globalShifter*site_index + catShifer*catCount,
                * vecPointer         = currentStateVector;
        
        for (long nodeCount = 0L; nodeCount<flatCLeaves.lLength; nodeCount++) {
            hyFloat *leafVec     = ((_CalcNode*)(((BaseRef*)flatCLeaves.list_data)[nodeCount]))->theProbs;
            CopyArray(vecPointer, leafVec, cBase);
            vecPointer += cBase;
        }
        
        // TODO SLKP 20180703: ugly fix for underflow which WON'T work if category count > 1
        
        for (long iNodeCount = 0L; iNodeCount < flatTree.lLength - 1; iNodeCount++) {
            node<long>* thisINode       = (node<long>*)flatNodes.list_data[iNodeCount];
            
            hyFloat sum = 0.;
            
            for (long cc = 0L; cc < cBase; cc++) {
                hyFloat      tmp = 1.0;
                
                for (long nc = 0; nc < thisINode->get_num_nodes(); nc++) {
                    hyFloat  tmp2 = 0.0;
                    _CalcNode   * child         = map_node_to_calcnode(thisINode->go_down(nc+1));
                    
                    hyFloat     * childSupport  = currentStateVector + node_to_index.GetDataByKey(child) *cBase,
                                * transMatrix   = child->GetCompExp(categoryCount>1?catCount:(-1))->theData + cc*cBase;
                    
                    for (long cc2 = 0; cc2 < cBase; cc2++) {
                        tmp2 += transMatrix[cc2] * childSupport[cc2];
                    }
                    
                    tmp *= tmp2;
                }
                vecPointer[cc] = tmp;
                sum += tmp;
            }
            
            if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                for (long cc = 0; cc < cBase; cc++) {
                    vecPointer[cc] *= _lfScalerUpwards;
                }
            }
            
            vecPointer += cBase;
            
        }
        RecoverNodeSupportStates2 (&GetRoot(),currentStateVector+globalShifter,currentStateVector,categoryCount>1?catCount:(-1), node_to_index);
    }
    /* pass 2; populate bottom-up vectors */
    /* for this we need to traverse the tree pre-order */
    /* because speed is not much of a concern, use a recursive call for compactness */
    
}

//_______________________________________________________________________________________________

void     _TheTree::RecoverNodeSupportStates2 (node<long>* thisNode, hyFloat * resultVector, hyFloat* forwardVector, long catID, _AVLListX& lookup) {
    
    _CalcNode   * thisNodeC     = map_node_to_calcnode (thisNode);
    hyFloat     * vecPointer    = resultVector + lookup.GetDataByKey(thisNodeC) * cBase;
    
    if (thisNode->parent) {
        if (thisNode->parent->parent) {
            hyFloat sum = 0.;
            for (long cc = 0; cc < cBase; cc++) {
                hyFloat tmp = 1.0;
                for (long nc = 0; nc < thisNode->parent->get_num_nodes(); nc++) {
                    hyFloat  tmp2            = 0.0;
                    _CalcNode   * child         = map_node_to_calcnode(thisNode->parent->go_down (nc+1));
                    bool          invert        = (child == thisNodeC);;
                    if (invert) {
                        child = map_node_to_calcnode (thisNode->parent);
                    }
                    
                    hyFloat  * childSupport  = invert?resultVector + cBase*lookup.GetDataByKey(child)
                    :forwardVector + lookup.GetDataByKey(child)*cBase,
                    * transMatrix   = child->GetCompExp(catID)->theData + cc*cBase;
                    
                    for (long cc2 = 0; cc2 < cBase; cc2++) {
                        tmp2 += transMatrix[cc2] * childSupport[cc2];
                    }
                    
                    tmp *= tmp2;
                }
                vecPointer[cc] = tmp;
                sum += tmp;
            }
            if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                for (long cc = 0; cc < cBase; cc++) {
                    vecPointer[cc] *= _lfScalerUpwards;
                }
            }
            
        } else {
            for (long cc = 0; cc < cBase; cc++,vecPointer++) {
                hyFloat tmp = 1.0;
                for (long nc = 0; nc < thisNode->parent->get_num_nodes(); nc++) {
                    hyFloat  tmp2            = 0.0;
                    _CalcNode   * child         = ((_CalcNode*)((BaseRef*)variablePtrs.list_data)[thisNode->parent->get_node(nc+1)->in_object]);
                    if (child != thisNodeC) {
                        hyFloat  * childSupport  = forwardVector + lookup.GetDataByKey(child)*cBase,
                        * transMatrix   = child->GetCompExp(catID)->theData + cc*cBase;
                        
                        for (long cc2 = 0; cc2 < cBase; cc2++) {
                            tmp2 += transMatrix[cc2] * childSupport[cc2];
                        }
                        
                        tmp *= tmp2;
                    }
                }
                *vecPointer = tmp;
            }
        }
    } else {
        InitializeArray (vecPointer, cBase, 1.0);
    }
    
    for (long nc = 0; nc < thisNode->get_num_nodes(); nc++) {
        RecoverNodeSupportStates2 (thisNode->get_node(nc+1),resultVector,forwardVector,catID,lookup);
    }
}


//_______________________________________________________________________________________________
_AVLListX*  _TheTree::ConstructNodeToIndexMap (bool doINodes) const {
    _SimpleList * nodes  = new _SimpleList;
    const _SimpleList * whichL = doINodes?&flatNodes:&flatLeaves;
    _AVLListX   * result = new _AVLListX (nodes);

    for (unsigned long   pistolero = 0UL; pistolero < whichL->lLength; pistolero++) {
        result->Insert ((BaseRef)whichL->list_data[pistolero], pistolero, false);
    }

    return        result;

}

  //_______________________________________________________________________________________________
void _TheTree::MapPostOrderToInOrderTraversal (_SimpleList& storeHere, bool doINodes) const {
  _AVLListX*          nodeMapper    = ConstructNodeToIndexMap (doINodes);

  _TreeIterator       ti (this, doINodes ? _HY_TREE_TRAVERSAL_PREORDER : _HY_TREE_TRAVERSAL_POSTORDER);

  unsigned long                allNodeCount = 0UL;

  storeHere.Populate (doINodes?flatTree.lLength:flatLeaves.lLength, 0, 0);

  while (_CalcNode* iterator = ti.Next()) {
    bool isTip = ti.IsAtLeaf();
    if ( isTip && !doINodes  || !isTip && doINodes) {
      storeHere.list_data[nodeMapper->GetXtra (nodeMapper->Find((BaseRef)(ti.GetNode())))] = allNodeCount++;
    }
  }

  nodeMapper->DeleteAll(false);
  DeleteObject (nodeMapper);
}

  //_______________________________________________________________________________________________
void    _TheTree::MarkDone (void) {
  _TreeIterator       ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  while (_CalcNode* iterator = ti.Next()) {
    iterator -> MarkDone();
  }
}

//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCostChar (_DataSetFilter const* dsf, long firstIndex, long secondIndex, const _SimpleList* child_count) const {

    const char *pastState = dsf->GetColumn(firstIndex),
               *thisState = dsf->GetColumn(secondIndex);

    long rootIndex = flatTree.lLength - 1,
         theCost = child_count->list_data[rootIndex],
         offset = flatLeaves.countitems();

    bool * marked_nodes = (bool*) calloc ( flatTree.lLength, sizeof(bool));
    //InitializeArray(marked_nodes, flatTree.lLength, false);
    //memset (marked_nodes, 0, sizeof(bool) * flatTree.lLength);
    
    for (long node_index = 0L; node_index < flatLeaves.lLength; node_index++) {
         long seq_index = dsf->theNodeMap.list_data[node_index];
         if (thisState[seq_index] != pastState[seq_index]) {
             /*long parentIndex = flatParents.list_data[node_index];
             while (marked_nodes[parentIndex] == false && parentIndex < rootIndex) {
                 marked_nodes[parentIndex]  = true;
                 theCost += child_count->list_data [parentIndex];
                 parentIndex = flatParents.list_data[parentIndex+offset];
             }*/
             marked_nodes [flatParents.list_data[node_index]] = true;
         }
    }

    // the root will always be changed, so don't need to check it
    for (long i=0; i<flatTree.lLength - 1; i++) {
        if (marked_nodes[i]) {
            marked_nodes [flatParents.list_data [i + offset]] = true;
            theCost += child_count->list_data [i];
        }
    }
    free (marked_nodes);

    return theCost;

}

//_______________________________________________________________________________________________

void    _TheTree::ClearConstraints (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        iterator->ClearConstraints();
    }
}


//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCost (_DataSetFilter const* dsf, long firstIndex, long secondIndex, _SimpleList* traversalTags, long orderIndex, _SimpleList const* childCount) const {

    long        filterL = dsf->GetPatternCount();

    
    
    bool * marked_nodes = (bool*) calloc (flatTree.lLength, sizeof(bool));
    //memset (marked_nodes, 0, sizeof(bool) * flatTree.lLength);
    //InitializeArray(marked_nodes, flatTree.lLength, false);

    
    dsf->CompareTwoSitesCallback (firstIndex, secondIndex, [marked_nodes, this] (long idx, unsigned long di)->void{
        marked_nodes [this->flatParents.list_data[di]] = true;
    });
    
    long theCost = 0,
         rootIndex = flatTree.lLength - 1;

    // don't check the root; it is always "dirty"
    
    if (childCount) {
        if (traversalTags && orderIndex) {
            for (long i=0; i<rootIndex; i++) {
               if (marked_nodes[i]) {
                   marked_nodes[flatParents.list_data[flatLeaves.lLength + i]] = true;
                   theCost += childCount->list_data[i];
               } else {
                   long theIndex = filterL * i + orderIndex;
                   traversalTags->list_data[theIndex/_HY_BITMASK_WIDTH_] |= bitMaskArray.masks[theIndex%_HY_BITMASK_WIDTH_];
               }
           }
        }
        
        for (long i=0; i<rootIndex; i++) {
            if (marked_nodes[i]) {
                marked_nodes[flatParents.list_data[flatLeaves.lLength + i]] = true;
                theCost += childCount->list_data[i];
            }
        }
    } else {
        for (long i=0; i<rootIndex; i++) {
            if (marked_nodes[i]) {
                marked_nodes[flatParents.list_data[flatLeaves.lLength + i]] = true;
                theCost += ((node <long>*)(flatNodes.list_data[i]))->get_num_nodes();
            } else if (traversalTags && orderIndex) {
                long theIndex = filterL * i + orderIndex;
                traversalTags->list_data[theIndex/_HY_BITMASK_WIDTH_] |= bitMaskArray.masks[theIndex%_HY_BITMASK_WIDTH_];
            }
        }
    }

    theCost += ((node <long>*)(flatNodes.list_data[rootIndex]))->get_num_nodes();
    
    free (marked_nodes);
    
    return theCost;

}


//_______________________________________________________________________________________________

void    _TheTree::MolecularClock (_String const& baseNode, _List& varsToConstrain) const {
    // TODO SLKP 20180803 : this needs a  review
    
    node<long>* topNode = nil;

    if (baseNode.empty()) { // called Molecular Clock on the entire tree
        topNode = &GetRoot();
        _String*  childNameP;
        if (rooted == ROOTED_LEFT) { // run separate constraint on the right child of the root
           MolecularClock (map_node_to_calcnode(theRoot->go_down (theRoot->get_num_nodes()))->ContextFreeName(), varsToConstrain);

        } else if (rooted == ROOTED_RIGHT) {
           MolecularClock (map_node_to_calcnode(theRoot->go_down (1))->ContextFreeName(), varsToConstrain);
        }
    } else {
        topNode = FindNodeByName (&baseNode);
    }

    if (!topNode) {
        HandleApplicationError (_String ("Molecular clock constraint has failed, since node '")
                   &baseNode
                   &"' is not a part of tree '"
                   &*GetName() & "'");
    } else
        for (unsigned long k=1UL; k<varsToConstrain.lLength; k++) {
            long varIndex = LocateVarByName (*(_String*)varsToConstrain (k));

            if (varIndex<0) {
                HandleApplicationError (_String ("Molecular clock constraint has failed, since variable' ") &*(_String*)varsToConstrain (k) &"' is undefined.");
                return ;
            }
            map_node_to_calcnode(topNode)->RecurseMC (variableNames.GetXtra(varIndex), topNode, true, rooted);
        }
}


  //_______________________________________________________________________________________________

void     _TheTree::AddNodeNamesToDS (_DataSet* ds, bool doTips, bool doInternals, char dOrS) const {
    // TODO SLKP 20180803 needs review
  if (dOrS == 2 && doTips && doInternals) {
    AddNodeNamesToDS (ds, false, true, 0);
    AddNodeNamesToDS (ds, true, false, 0);
    return;
  }

  RetrieveNodeNames (doTips, doInternals,dOrS ? _HY_TREE_TRAVERSAL_POSTORDER : _HY_TREE_TRAVERSAL_PREORDER).ForEach ([&ds] (BaseRef name, unsigned long )->void {
        ds->AddName (*(_String*)name);
    });
}

// LF COMPUTE FUNCTIONS
// TODO SLKP 20180803 these all could use a review

//extern long mes_counter;

/*----------------------------------------------------------------------------------------------------------*/
void        _TheTree::ExponentiateMatrices  (_List& expNodes, long tc, long catID) {
    _List           matrixQueue, nodesToDo;
    
    _SimpleList     isExplicitForm ((unsigned long)expNodes.countitems());
    bool            hasExpForm = false;
    
    //long b4 = mes_counter;

    _SimpleList dpv;
    _keepTrackOfDepVars = new _AVLList (&dpv);
    useGlobalUpdateFlag = true;
    
    for (unsigned long nodeID = 0; nodeID < expNodes.lLength; nodeID++) {
        long didIncrease = matrixQueue.lLength;
        _CalcNode* thisNode = (_CalcNode*) expNodes(nodeID);
        //if (hasExpForm) {
        //    b4 = mes_counter;
        //}
        if (thisNode->RecomputeMatrix (catID, categoryCount, nil, &matrixQueue,&isExplicitForm)) {
            hasExpForm = true;
        }
        //if (hasExpForm && mes_counter - b4 > 3) {
        //    printf ("%d/%d %s\n", mes_counter - b4, likeFuncEvalCallCount , thisNode->GetName()->get_str());
        //}
#ifdef _UBER_VERBOSE_DUMP
        if (likeFuncEvalCallCount == _UBER_VERBOSE_DUMP)
            printf ("NodeID %ld (%s). Old length %ld, new length %ld (%ld)\n", nodeID, thisNode->GetName()->sData, didIncrease,matrixQueue.lLength, isExplicitForm.lLength);
#endif
        if ((didIncrease = (matrixQueue.lLength - didIncrease))) {
            for (long copies = 0; copies < didIncrease; copies++) {
                nodesToDo << thisNode;
            }
        }
    }
    
    //ObjectToConsole(_keepTrackOfDepVars);
    //NLToConsole();
    ResetDepComputedFlags (dpv);
    DeleteAndZeroObject(_keepTrackOfDepVars);
    
    useGlobalUpdateFlag = false;

        
    _List * computedExponentials = hasExpForm? new _List (matrixQueue.lLength) : nil;
    
    _SimpleList parallel, serial;
    isExplicitForm.Each ([&parallel,&serial](long mx_count, unsigned long id) -> void {
        if (mx_count < 0) serial << id; else parallel << id;
    });
    hy_global::matrix_exp_count += matrixQueue.lLength - serial.countitems();
    unsigned long nt = cBase<20 ? 1:(MIN(tc, (cBase < 60 ? parallel.lLength >> 1 : parallel.lLength) + 1));
    unsigned long cs = cBase<20 ? 10 : (cBase < 60 ? 2 : 1);

    //printf ("_TheTree::ExponentiateMatrices %d total, %d no update, (block update %d)\n", parallel.lLength, serial.lLength, nt);

    if (parallel.lLength) {
#ifdef _OPENMP
  #if _OPENMP>=201511
    #pragma omp parallel for default(shared) schedule(dynamic, cs) proc_bind(spread) if (nt>1)  num_threads (nt)
  #else
  #if _OPENMP>=200803
    #pragma omp parallel for default(shared) schedule(guided) proc_bind(spread) if (nt>1)  num_threads (nt)
  #endif
#endif
#endif
        for  (long id = 0L; id < parallel.lLength; id++) {
            long matrixID = parallel.get (id);
            if (isExplicitForm.list_data[matrixID] == 0 || !hasExpForm) { // normal matrix to exponentiate
                ((_CalcNode*) nodesToDo(matrixID))->SetCompExp ((_Matrix*)matrixQueue(matrixID), catID, true);
            } else {
                (*computedExponentials) [matrixID] = ((_Matrix*)matrixQueue(matrixID))->Exponentiate(1., true);
            }
        }
    }
    
    for (long id = 0L; id < serial.lLength; id++) {
        long matrixID = serial.get (id);
        _Matrix *already_computed = ((_Matrix*)matrixQueue(matrixID));
        (*computedExponentials) [matrixID] = already_computed;
    }
    
 
    if (computedExponentials) {
        _CalcNode * current_node         = nil;
        _List       buffered_exponentials;
        
        for (unsigned long mx_index = 0; mx_index < nodesToDo.lLength; mx_index++) {
            if (isExplicitForm.list_data[mx_index] > 0) {
                _CalcNode *next_node = (_CalcNode*) nodesToDo (mx_index);
                 if (next_node != current_node) {
                    if (current_node) {
                        current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
                    }
                    current_node = next_node;
                    buffered_exponentials.Clear();
                    buffered_exponentials << (*computedExponentials)(mx_index);
                }
                else {
                    buffered_exponentials << (*computedExponentials)(mx_index);
                }
            } else {
                if (isExplicitForm.list_data[mx_index] < 0) {
                    if (current_node && buffered_exponentials.countitems()) {
                        current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
                    }
                    buffered_exponentials.Clear();
                    buffered_exponentials << (*computedExponentials)(mx_index);
                    ((_CalcNode*) nodesToDo (mx_index))->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials, true);
                    buffered_exponentials.Clear();
                } else {
                    if (current_node) {
                        current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
                    }
                }
                current_node = nil;
            }
        }
        if (current_node) {
            current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
        }
        DeleteObject(computedExponentials);
#ifdef _UBER_VERBOSE_DUMP_MATRICES
        if (likeFuncEvalCallCount == _UBER_VERBOSE_DUMP) {
            fprintf (stderr, "\n T_MATRIX = {");
            for (unsigned long nodeID = 0; nodeID < flatLeaves.lLength + flatTree.lLength - 1; nodeID++) {
                bool    isLeaf     = nodeID < flatLeaves.lLength;
                
                _CalcNode * current_node = isLeaf? (((_CalcNode**) flatCLeaves.list_data)[nodeID]):
                (((_CalcNode**) flatTree.list_data)  [nodeID - flatLeaves.lLength]);
                if (nodeID) {
                    fprintf (stderr, ",");
                }
                fprintf (stderr, "\n\"%s\":%s", current_node->GetName()->get_str(), _String((_String*)current_node->GetCompExp()->toStr()).get_str());
                
            }
            fprintf (stderr, "\n};\n");
        }
#endif
        
        
    }
}

/*----------------------------------------------------------------------------------------------------------*/

long        _TheTree::DetermineNodesForUpdate   (_SimpleList& updateNodes, _List* expNodes, long catID, long addOne, bool canClear,
                                                _AVLListX * var_to_node_map, _AVLList * changed_variables) {
  _CalcNode       *currentTreeNode;
  long            lastNodeID = -1;
  unsigned long   tagged_node_count = 0UL;
    
    _SimpleList     _ctSupp;
    _AVLListX       ctCache (&_ctSupp);

    // look for nodes with model changes and mark the path up to the root as needing an update
        
  #define DIRECT_INDEX(N) (flatParents.list_data[N]+flatLeaves.lLength)
    
  auto _handle_node = [&] (long node_id, bool do_list)->void {
      bool    isLeaf     = node_id < flatLeaves.lLength;
      
      if (isLeaf) {
        currentTreeNode = (((_CalcNode**) flatCLeaves.list_data)[node_id]);
      } else {
        currentTreeNode = (((_CalcNode**) flatTree.list_data)  [node_id - flatLeaves.lLength]);
      }
      
      if (currentTreeNode->NeedNewCategoryExponential (catID, &ctCache)) {
        if (expNodes) {
          (*expNodes) << currentTreeNode;
          lastNodeID = node_id;
        } else {
          currentTreeNode->RecomputeMatrix (catID, categoryCount, nil);
        }
        
        if (do_list) {
            nodesToUpdate.list_data[node_id] = 2;
            tagged_node_count++;
        }
      }
      
      if (do_list && nodesToUpdate.list_data[node_id]) {
        long parent_index = DIRECT_INDEX(node_id);
        //if (nodesToUpdate.list_data[parent_index] == 0) {
            nodesToUpdate.list_data[parent_index] = 2;
            tagged_node_count++;
        //}
      }
  };

  bool already_done = false;
  if (changed_variables && var_to_node_map && changed_variables->countitems() < 8 && forceRecalculationOnTheseBranches.empty()) {
      already_done = true;
      //_SimpleList   nodes;
      updateNodes.RequestSpace (changed_variables->countitems());
      _AVLList uniques (&updateNodes);

      for (long i = 0L; i < changed_variables->dataList->lLength; i++) {
          long var_index = changed_variables->dataList->list_data[i];
          long node_index = var_to_node_map->FindAndGetXtra ((BaseRefConst)var_index,-1L);
          if (node_index < 0) {
              already_done = false;
              break;
          }
          //updateNodes << node_index;
          uniques.InsertNumber(node_index);
      }
      //already_done = false;
      //if (likeFuncEvalCallCount == 7) {
      //    ObjectToConsole(&updateNodes); NLToConsole();
      //}
      if (already_done) {
          if (addOne >= 0) {
              uniques.InsertNumber(addOne);
          }
          long i = 0L;
          for (; i < updateNodes.lLength; i++) {
              _handle_node (updateNodes.get (i), false);
          }
          //if (likeFuncEvalCallCount == 7) {
          //    ObjectToConsole(&updateNodes); NLToConsole();
          //}
          

          for (long i2 = 0; i2 < i; i2++) {
              long  my_index = updateNodes.get (i2),
                    parent_index = flatParents.list_data[my_index];
                            
              while (parent_index >= 0) {
                  long dir_index = parent_index + flatLeaves.lLength;
                  //if (uniques.InsertNumber(dir_index) >= 0) { // performed insertion
                  //    node<long>* parent_node = (node<long>*)flatNodes.get (parent_index);
                  //}
                  uniques.InsertNumber(dir_index);
                  parent_index = flatParents.list_data[dir_index];
              }
          }

          _SimpleList children;
          
          /*nodesToUpdate.Populate (flatTree.lLength, 0, 0);
          
          for (long i = 0; i < updateNodes.lLength; i++) {
              long tagged_node = updateNodes.get (i);
              if (tagged_node >= flatLeaves.lLength) {
                  nodesToUpdate.list_data [tagged_node - flatLeaves.lLength] = true;
              }
          }
          
          for (long i = 0; i < flatParents.lLength - 1; i++) {
              long my_parent = flatParents.list_data [i];
              if (nodesToUpdate.list_data[my_parent]) {
                  children << i;
              }
          }*/
          
          for (long i = 0; i < flatParents.lLength - 1; i++) {
              long my_parent = flatParents.list_data [i] + flatLeaves.lLength;
              if (uniques.FindLong (my_parent) >= 0) {
                  children << i;
              }
          }
          
          for (long i = 0; i < children.lLength; i++) {
              uniques.InsertNumber(children.get (i));
          }
          
           if (uniques.countitems() <= 32) {
              updateNodes.BubbleSort();
          } else {
              uniques.ReorderList();
          }
          
          updateNodes.Pop();
          
      } else {
          updateNodes.Clear(false);
      }
  }
    
    //bool frc = forceRecalculationOnTheseBranches.nonempty();

    if (!already_done) {
       nodesToUpdate.Populate (flatLeaves.lLength + flatTree.lLength, 0, 0);
       
       if (addOne >= 0) {
         nodesToUpdate.list_data[addOne] = 2;
         tagged_node_count++;
       }
       
       if (forceRecalculationOnTheseBranches.nonempty()) {
         forceRecalculationOnTheseBranches.Each ([this] (long value, unsigned long) -> void {
           this->nodesToUpdate.list_data [value] = 2L;
         });
         tagged_node_count += forceRecalculationOnTheseBranches.countitems();
         
         if (canClear) {
           forceRecalculationOnTheseBranches.Clear();
         }
       }

       //if (frc) {        ObjectToConsole(&nodesToUpdate); NLToConsole();}
       for (unsigned long nodeID = 0UL; nodeID < nodesToUpdate.lLength - 1UL; nodeID++) {
            _handle_node (nodeID, true);
       }
       //if (frc) {        ObjectToConsole(&nodesToUpdate); NLToConsole();}

      
        // one more pass to pick up all DIRECT descendants of changed internal nodes
      
       for (unsigned long nodeID = 0UL; nodeID < nodesToUpdate.lLength - 1UL; nodeID++)
        if (nodesToUpdate.list_data[nodeID] == 0 && nodesToUpdate.list_data[DIRECT_INDEX(nodeID)] == 2) {
          nodesToUpdate.list_data[nodeID] = 1;
          tagged_node_count++;
        }
       
       //if (frc) {        ObjectToConsole(&nodesToUpdate); NLToConsole();}

        // write out all changed nodes
      updateNodes.RequestSpace(tagged_node_count);
      // 20200610: for larger trees, this helps with reducing memory allocations
      for (unsigned long nodeID = 0UL; nodeID < nodesToUpdate.lLength - 1UL; nodeID++) {
        if (nodesToUpdate.list_data[nodeID]) {
          updateNodes << nodeID;
            /*if (frc) {
                _CalcNode       *currentTreeNode;
                if (nodeID < flatLeaves.lLength) {
                    currentTreeNode = ((_CalcNode*) flatCLeaves (nodeID));
                } else {
                    currentTreeNode = ((_CalcNode*) flatTree    (nodeID-flatLeaves.lLength));
                }
                printf ("%d %s\n", nodeID, currentTreeNode->GetName()->get_str());
            }*/
        }
      }
   }

    

   /*printf ("%ld ", likeFuncEvalCallCount);
   ObjectToConsole(&updateNodes);
   printf ("\n");*/

  if (expNodes && expNodes->countitems() == 1) {
    return lastNodeID;
  }
  
  return -1;
}

/*----------------------------------------------------------------------------------------------------------*/

void        _TheTree::FillInConditionals        (_DataSetFilter const*        theFilter, hyFloat*  iNodeCache,  _SimpleList*   tcc)
// this utility function will simply fill in all the conditional probability vectors for internal nodes,
// including those that were skipped due to column sorting optimization
// this is useful to avoid code duplication for other functions (e.g. ancestral sampling) that
// make use of conditional probability vectors, but would not benefit from subtree caching
{
    if (!tcc) {
        return;
    }
    
    long            alphabetDimension     =         theFilter->GetDimension(),
    siteCount           =         theFilter->GetPatternCount();
    
    for  (long nodeID = 0; nodeID < flatTree.lLength; nodeID++) {
        hyFloat * conditionals       = iNodeCache +(nodeID  * siteCount) * alphabetDimension;
        long        currentTCCIndex     = siteCount * nodeID,
        currentTCCBit        = currentTCCIndex % _HY_BITMASK_WIDTH_;
        
        currentTCCIndex /= _HY_BITMASK_WIDTH_;
        for (long siteID = 0; siteID < siteCount; siteID++, conditionals += alphabetDimension) {
            if (siteID  && (tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0) {
                for (long k = 0; k < alphabetDimension; k++) {
                    conditionals[k] = conditionals[k-alphabetDimension];
                }
            }
            if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }
            
        }
    }
}




/*----------------------------------------------------------------------------------------------------------*/

const _CalcNode* _TheTree::GetNodeFromFlatIndex(long index) const {
    return index < flatLeaves.lLength ? (((_CalcNode**) flatCLeaves.list_data)[index]):
    (((_CalcNode**) flatTree.list_data)   [index - flatLeaves.lLength]);
}

/*----------------------------------------------------------------------------------------------------------*/

hyFloat          _TheTree::ComputeLLWithBranchCache (
                                                     _SimpleList&            siteOrdering,
                                                     long                    brID,
                                                     hyFloat*         cache,
                                                     _DataSetFilter const*     theFilter,
                                                     long                    siteFrom,
                                                     long                    siteTo,
                                                     long                    catID,
                                                     hyFloat*         storageVec
                                                     )
{
    auto bookkeeping =  [&siteOrdering, &storageVec, &theFilter] (const long siteID, const hyFloat accumulator, hyFloat& correction, hyFloat& result) ->  void {
        
        long direct_index = siteOrdering.list_data[siteID];
        
        if (storageVec) {
            storageVec [direct_index] = accumulator;
        } else {
            if (accumulator <= 0.0) {
                
                
                /*fprintf (stderr, "ZERO TERM AT SITE %ld (direct %ld) EVAL %ld\n",siteID,direct_index, likeFuncEvalCallCount);
                for (long s = 0; s < theFilter->NumberSpecies(); s++) {
                    fprintf (stderr, "%s:%s\n", theFilter->GetSequenceName(s)->get_str(), theFilter->RetrieveState(direct_index, s).get_str());
                }
                fprintf (stderr, "\n");*/
                
                throw (1L+direct_index);
            }
            
            hyFloat term;
            long       site_frequency = theFilter->theFrequencies.get(direct_index);
            if ( site_frequency > 1L) {
                term =  log(accumulator) * site_frequency - correction;
            } else {
                term = log(accumulator) - correction;
            }
            
            
            /*if (likeFuncEvalCallCount == 328 && siteID == 30) {
                fprintf (stderr, "CACHE, %ld, %ld, %20.15lg, %20.15lg, %20.15lg,  %20.15lg\n", likeFuncEvalCallCount, siteID, accumulator, correction, term, result);
            }*/
            
            
            hyFloat temp_sum = result + term;
            correction = (temp_sum - result) - term;
            result = temp_sum;
            //result += log(accumulator) * theFilter->theFrequencies [siteOrdering.list_data[siteID]];
        }
    };
    
    const unsigned long          alphabetDimension      = theFilter->GetDimension(),
    siteCount              =  theFilter->GetPatternCount();
    
    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }
    
    hyFloat const * branchConditionals = cache              + siteFrom * alphabetDimension;
    hyFloat const * rootConditionals   = branchConditionals + siteCount * alphabetDimension;
    hyFloat  result = 0.0,
    correction = 0.0;
    
    
    //printf ("ComputeLLWithBranchCache @ %d catID = %d branchID = %d\n", likeFuncEvalCallCount, catID, brID);
    
    _CalcNode const *givenTreeNode = GetNodeFromFlatIndex (brID);
    
    hyFloat  const * transitionMatrix = givenTreeNode->GetCompExp(catID)->theData;
    
    
    // cases by alphabet dimension
    /*if (likeFuncEvalCallCount == 15098) {
        fprintf (stderr, "\nBRANCH ID %ld (%s)\n", brID, givenTreeNode->GetName()->get_str());
        for (long e = 0; e < 16; e++) {
            fprintf (stderr, "%ld => %lg\n", transitionMatrix[e]);
        }
    }*/
    
    try {
        switch (alphabetDimension) {
                /****
                 
                 NUCLEOTIDES
                 
                 ****/
            case 4UL: {
#ifdef _SLKP_USE_AVX_INTRINSICS
                __m256d tmatrix_transpose [4] = {
                    (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                    (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                    (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                    (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                };
#endif
                
#ifdef _SLKP_USE_ARM_NEON
                float64x2x2_t tmatrix_transpose [4] = {
                    (float64x2x2_t) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                    (float64x2x2_t) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                    (float64x2x2_t) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                    (float64x2x2_t) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                };
#endif
                
                for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
                    hyFloat accumulator = 0.;
                    /*if (likeFuncEvalCallCount == 328 && siteID == 30) {
                        fprintf (stderr, "CACHE, %ld, %ld, %20.15lg/%20.15lg, %20.15lg/%20.15lg, %20.15lg/%20.15lg, %20.15lg/%20.15lg\n", likeFuncEvalCallCount, siteID, rootConditionals[0],branchConditionals[0], rootConditionals[1], branchConditionals[1], rootConditionals[2], branchConditionals[2], rootConditionals[3], branchConditionals[3]);
                        
                        fprintf (stderr, "\n%s (%ld)\n", givenTreeNode->GetName()->get_str(), brID);
                        for (int i = 0; i < 16; i++) {
                            fprintf (stderr, "%ld\t%20.15lg\n", i, transitionMatrix[i]);
                        }
                    }*/

#ifdef _SLKP_USE_AVX_INTRINSICS

        
                    __m256d     c0     = _mm256_set1_pd  (branchConditionals[0]-branchConditionals[3]),
                                c1     = _mm256_set1_pd  (branchConditionals[1]-branchConditionals[3]),
                                c2     = _mm256_set1_pd  (branchConditionals[2]-branchConditionals[3]),
                                c3     = _mm256_set1_pd  (branchConditionals[3]);
    
                    __m256d     t[2];
  
    
                    t[0] = _hy_matrix_handle_axv_mfma (_mm256_mul_pd (c0, tmatrix_transpose[0]), c1, tmatrix_transpose[1]);
                    t[1] = _hy_matrix_handle_axv_mfma (c3,c2, tmatrix_transpose[2]);
    
                    t[0] = _mm256_mul_pd (_mm256_loadu_pd (rootConditionals), _mm256_add_pd (t[0], t[1]));
    
                    accumulator = _avx_sum_4 (_mm256_mul_pd (t[0], _mm256_loadu_pd (theProbs)));
        

                    /*__m256d root_c = _mm256_loadu_pd (rootConditionals),
                    probs  = _mm256_loadu_pd (theProbs),
                    b_cond0 = _mm256_set1_pd(branchConditionals[0]),
                    b_cond1 = _mm256_set1_pd(branchConditionals[1]),
                    b_cond2 = _mm256_set1_pd(branchConditionals[2]),
                    b_cond3 = _mm256_set1_pd(branchConditionals[3]),
                    s01    = _mm256_add_pd ( _mm256_mul_pd (b_cond0, tmatrix_transpose[0]), _mm256_mul_pd (b_cond1, tmatrix_transpose[1])),
                    s23    = _mm256_add_pd ( _mm256_mul_pd (b_cond2, tmatrix_transpose[2]), _mm256_mul_pd (b_cond3, tmatrix_transpose[3]));
                    accumulator = _avx_sum_4(_mm256_mul_pd (_mm256_mul_pd (root_c, probs), _mm256_add_pd (s01,s23)));*/
                    
#elif defined _SLKP_USE_ARM_NEON
                
                    float64x2_t c0     = vdupq_n_f64(branchConditionals[0]-branchConditionals[3]),
                                c1     = vdupq_n_f64(branchConditionals[1]-branchConditionals[3]),
                                c2     = vdupq_n_f64(branchConditionals[2]-branchConditionals[3]),
                                c3     = vdupq_n_f64(branchConditionals[3]);
                    
                    float64x2_t t[4];
                  
                    t[0] = vfmaq_f64 (vmulq_f64  (c1,tmatrix_transpose[1].val[0]),c0,tmatrix_transpose[0].val[0]);
                    t[1] = vfmaq_f64 (vmulq_f64  (c1,tmatrix_transpose[1].val[1]),c0,tmatrix_transpose[0].val[1]);
                    t[2] = vaddq_f64 (vmulq_f64  (c2,tmatrix_transpose[2].val[0]),c3);
                    t[3] = vaddq_f64 (vmulq_f64  (c2,tmatrix_transpose[2].val[1]),c3);
                    
                    t[0] = vmulq_f64  (vld1q_f64 (rootConditionals), vaddq_f64 (t[0], t[2]));
                    t[1] = vmulq_f64  (vld1q_f64 (rootConditionals+2), vaddq_f64 (t[1],t[3]));
                    t[0] = vmulq_f64 (t[0], vld1q_f64 (theProbs));
                    t[1] = vmulq_f64 (t[1], vld1q_f64 (theProbs+2));
                    accumulator = _neon_sum_2 (vaddq_f64 (t[0],t[1]));
                    

                    
#else
                    /*accumulator =    rootConditionals[0] * theProbs[0] *
                    (branchConditionals[0] *  transitionMatrix[0] + branchConditionals[1] *  transitionMatrix[1] + branchConditionals[2] *  transitionMatrix[2] + branchConditionals[3] *  transitionMatrix[3]) +
                    rootConditionals[1] * theProbs[1] *
                    (branchConditionals[0] *  transitionMatrix[4] + branchConditionals[1] *  transitionMatrix[5] + branchConditionals[2] *  transitionMatrix[6] + branchConditionals[3] *  transitionMatrix[7]) +
                    rootConditionals[2] * theProbs[2] *
                    (branchConditionals[0] *  transitionMatrix[8] + branchConditionals[1] *  transitionMatrix[9] + branchConditionals[2] *  transitionMatrix[10] + branchConditionals[3] *  transitionMatrix[11]) +
                    rootConditionals[3] * theProbs[3] *
                    (branchConditionals[0] *  transitionMatrix[12] + branchConditionals[1] *  transitionMatrix[13] + branchConditionals[2] *  transitionMatrix[14] + branchConditionals[3] *  transitionMatrix[15]);
                    */
                    
                    hyFloat t1 = branchConditionals[0] - branchConditionals[3],
                    t2 = branchConditionals[1] - branchConditionals[3],
                    t3 = branchConditionals[2] - branchConditionals[3],
                    t4 = branchConditionals[3];
                    
                    accumulator =
                      rootConditionals[0] * theProbs[0] * ((transitionMatrix[0]  * t1 + transitionMatrix[1] * t2) + (transitionMatrix[2] * t3 + t4))
                    + rootConditionals[1] * theProbs[1] * ((transitionMatrix[4]  * t1 + transitionMatrix[5] * t2) + (transitionMatrix[6] * t3 + t4))
                    + rootConditionals[2] * theProbs[2] * ((transitionMatrix[8]  * t1 + transitionMatrix[9] * t2) + (transitionMatrix[10] * t3 + t4))
                    + rootConditionals[3] * theProbs[3] * ((transitionMatrix[12] * t1 + transitionMatrix[13] * t2) + (transitionMatrix[14] * t3 + t4));
                    
                    //return (parentConditionals[0] + parentConditionals[1]) + (parentConditionals[2] + parentConditionals[3]);
                    
                    
#endif
                    bookkeeping (siteID, accumulator, correction, result);
                    rootConditionals += 4UL;
                    branchConditionals += 4UL;
                 } // siteID
            }
                break;
                /****
                 
                 AMINOACIDS
                 
                 ****/
            case 20UL: {
                for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
                    hyFloat accumulator = 0.;
#ifdef _SLKP_USE_AVX_INTRINSICS
                    __m256d bc_vector[5] = {_mm256_loadu_pd(branchConditionals),
                        _mm256_loadu_pd(branchConditionals+4UL),
                        _mm256_loadu_pd(branchConditionals+8UL),
                        _mm256_loadu_pd(branchConditionals+12UL),
                        _mm256_loadu_pd(branchConditionals+16UL)};
                    
                    
                    hyFloat const * tm = transitionMatrix;
                    
                    for (unsigned long p = 0UL; p < 20UL; p++, rootConditionals++) {
                        
                        __m256d t_matrix[5] = { _mm256_loadu_pd(tm),
                            _mm256_loadu_pd(tm+4UL),
                            _mm256_loadu_pd(tm+8UL),
                            _mm256_loadu_pd(tm+12UL),
                            _mm256_loadu_pd(tm+16UL)};
                        
                        
                        t_matrix[0] = _mm256_mul_pd(t_matrix[0], bc_vector[0]);
                        t_matrix[1] = _mm256_mul_pd(t_matrix[1], bc_vector[1]);
                        t_matrix[2] = _mm256_mul_pd(t_matrix[2], bc_vector[2]);
                        t_matrix[3] = _mm256_mul_pd(t_matrix[3], bc_vector[3]);
                        t_matrix[4] = _mm256_mul_pd(t_matrix[4], bc_vector[4]);
                        
                        t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[1]);
                        t_matrix[1] = _mm256_add_pd (t_matrix[2],t_matrix[3]);
                        t_matrix[3] = _mm256_add_pd (t_matrix[0],t_matrix[1]);
                        
                        tm += 20UL;
                        
                        accumulator += *rootConditionals * theProbs[p] * _avx_sum_4(_mm256_add_pd (t_matrix[3],t_matrix[4]));
                    }
                    
#elif defined _SLKP_USE_ARM_NEON
                    
                    float64x2_t bc_vector[10] = {vld1q_f64(branchConditionals),
                        vld1q_f64(branchConditionals+2UL),
                        vld1q_f64(branchConditionals+4UL),
                        vld1q_f64(branchConditionals+6UL),
                        vld1q_f64(branchConditionals+8UL),
                        vld1q_f64(branchConditionals+10UL),
                        vld1q_f64(branchConditionals+12UL),
                        vld1q_f64(branchConditionals+14UL),
                        vld1q_f64(branchConditionals+16UL),
                        vld1q_f64(branchConditionals+18UL)
                    };
 
                    hyFloat const * tm = transitionMatrix;
                    
                    for (unsigned long p = 0UL; p < 20UL; p++, rootConditionals++) {
                        
                        float64x2_t S1 = vmulq_f64  ( bc_vector[0], vld1q_f64(tm)),
                                    S2 = vmulq_f64 ( bc_vector[1], vld1q_f64(tm+2UL));
                        
                        S1 = vfmaq_f64 (S1, bc_vector[2], vld1q_f64(tm+4UL));
                        S2 = vfmaq_f64 (S2, bc_vector[3], vld1q_f64(tm+6UL));
                        S1 = vfmaq_f64 (S1, bc_vector[4], vld1q_f64(tm+8UL));
                        S2 = vfmaq_f64 (S2, bc_vector[5], vld1q_f64(tm+10UL));
                        S1 = vfmaq_f64 (S1, bc_vector[6], vld1q_f64(tm+12UL));
                        S2 = vfmaq_f64 (S2, bc_vector[7], vld1q_f64(tm+14UL));
                        S1 = vfmaq_f64 (S1, bc_vector[8], vld1q_f64(tm+16UL));
                        S2 = vfmaq_f64 (S2, bc_vector[9], vld1q_f64(tm+18UL));
                        
                        accumulator += *rootConditionals * theProbs[p] * _neon_sum_2(vaddq_f64 (S1,S2));
                        tm += 20UL;
                   }
#else // _SLKP_USE_AVX_INTRINSICS
                    unsigned long rmx = 0UL;
                    for (unsigned long p = 0UL; p < 20UL; p++,rootConditionals++) {
                        hyFloat     r2 = 0.;
                        
                        for (unsigned long c = 0UL; c < 20UL; c+=4UL, rmx +=4UL) {
                            r2 += (branchConditionals[c]   *  transitionMatrix[rmx] +
                                   branchConditionals[c+1] *  transitionMatrix[rmx+1]) +
                            (branchConditionals[c+2] *  transitionMatrix[rmx+2] +
                             branchConditionals[c+3] *  transitionMatrix[rmx+3]);
                        }
                        
                        accumulator += *rootConditionals * theProbs[p] * r2;
                    }
#endif  // _SLKP_USE_AVX_INTRINSICS
                    branchConditionals += 20UL;
                    bookkeeping (siteID, accumulator, correction, result);
                    
                } // siteID
                
            } // case 20
                break;
                /****
                 
                 CODONS
                 
                 ****/
            case 60UL:
            case 61UL:
            case 62UL:
            case 63UL: {
                for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
                    hyFloat accumulator = 0.;
                    
                    unsigned long rmx = 0UL;
                    for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
                        hyFloat     r2 = 0.;
                        unsigned       long c = 0UL;
                        
#ifdef _SLKP_USE_AVX_INTRINSICS
                        
                        __m256d sum256 = _mm256_setzero_pd ();
                        
                        for (; c < 60UL; c+=12UL, rmx +=12UL) {
                            
                            __m256d branches0 = _mm256_loadu_pd (branchConditionals+c),
                            branches1 = _mm256_loadu_pd (branchConditionals+c+4),
                            branches2 = _mm256_loadu_pd (branchConditionals+c+8),
                            matrix0   = _mm256_loadu_pd (transitionMatrix+rmx),
                            matrix1   = _mm256_loadu_pd (transitionMatrix+rmx+4),
                            matrix2   = _mm256_loadu_pd (transitionMatrix+rmx+8);
                            
                            branches0 = _mm256_mul_pd(branches0, matrix0);
                            branches1 = _mm256_mul_pd(branches1, matrix1);
                            branches2 = _mm256_mul_pd(branches2, matrix2);
                            
                            branches0 = _mm256_add_pd (branches0,branches2);
                            sum256 = _mm256_add_pd (branches0,_mm256_add_pd (sum256, branches1));
                        }
                        
                        r2 = _avx_sum_4(sum256);
#elif defined _SLKP_USE_ARM_NEON
                        
                        float64x2_t sum128 = vdupq_n_f64 (0.);
                        for (; c < 60UL; c+=6UL, rmx +=6UL) {
                            
                            float64x2_t branches0 = vld1q_f64 (branchConditionals+c),
                                        branches1 = vld1q_f64 (branchConditionals+c+2),
                                        branches2 = vld1q_f64 (branchConditionals+c+4),
                                        matrix0   = vld1q_f64 (transitionMatrix+rmx),
                                        matrix1   = vld1q_f64 (transitionMatrix+rmx+2),
                                        matrix2   = vld1q_f64 (transitionMatrix+rmx+4);
                            
                            branches0 = vfmaq_f64 (vmulq_f64 (branches0, matrix0), branches1, matrix1);
                            branches1 = vfmaq_f64 (sum128, branches2, matrix2);
                            sum128 = vaddq_f64 (branches1, branches0);
                        }
                        r2 = _neon_sum_2(sum128);
                        
#else // _SLKP_USE_AVX_INTRINSICS
                        for (; c < 60UL; c+=4UL, rmx +=4UL) {
                            r2 += (branchConditionals[c]   *  transitionMatrix[rmx] +
                                   branchConditionals[c+1] *  transitionMatrix[rmx+1]) +
                            (branchConditionals[c+2] *  transitionMatrix[rmx+2] +
                             branchConditionals[c+3] *  transitionMatrix[rmx+3]);
                        }
#endif
                        
                        for (; c < alphabetDimension; c++, rmx ++) {
                            r2 += branchConditionals[c]   *  transitionMatrix[rmx];
                        }
                        
                        accumulator += *rootConditionals * theProbs[p] * r2;
                    }
                    
                    branchConditionals += alphabetDimension;
                    bookkeeping (siteID, accumulator, correction, result);
                    
                }
            } // cases 60-63
                break;
            default: { // valid alphabetDimension >= 2
                
                if (alphabetDimension % 2) { // odd
                    unsigned long alphabetDimension_minus1 = alphabetDimension-1;
                    for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
                        hyFloat accumulator = 0.;
                        
                        unsigned long rmx = 0UL;
                        for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
                            hyFloat     r2 = 0.;
                            
                            for (unsigned long c = 0UL; c < alphabetDimension_minus1; c+=2UL, rmx +=2UL) {
                                r2 +=  branchConditionals[c]   *  transitionMatrix[rmx] +
                                branchConditionals[c+1] *  transitionMatrix[rmx+1];
                            }
                            
                            r2 += branchConditionals[alphabetDimension_minus1]   *  transitionMatrix[rmx++];
                            
                            accumulator += *rootConditionals * theProbs[p] * r2;
                        }
                        
                        branchConditionals += alphabetDimension;
                        bookkeeping (siteID, accumulator, correction, result);
                        
                    }
                } else {
                    for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
                        hyFloat accumulator = 0.;
                        
                        unsigned long rmx = 0UL;
                        for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
                            hyFloat     r2 = 0.;
                            
                            for (unsigned long c = 0UL; c < alphabetDimension; c+=2UL, rmx +=2UL) {
                                r2 +=  branchConditionals[c]   *  transitionMatrix[rmx] +
                                branchConditionals[c+1] *  transitionMatrix[rmx+1];
                            }
                            
                            accumulator += *rootConditionals * theProbs[p] * r2;
                        }
                        
                        branchConditionals += alphabetDimension;
                        bookkeeping (siteID, accumulator, correction, result);
                        
                    }
                }
                
            } // default
        } // switch (alphabetDimension)
    } catch (long site) {
#pragma omp critical
        {
            hy_global::ReportWarning (_String("Site ") & _String(site) & " evaluated to a 0 probability in ComputeLLWithBranchCache");
        }
        return -INFINITY;
    }
    return result;
}

/*----------------------------------------------------------------------------------------------------------*/

hyFloat      _TheTree::ComputeTwoSequenceLikelihood
(
 _SimpleList   & siteOrdering,
 _DataSetFilter const* theFilter,
 long      *         lNodeFlags,
 _Vector* lNodeResolutions,
 long                siteFrom,
 long                siteTo,
 long                catID,
 hyFloat*     storageVec
 )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
    // process the leaves first
    
    long            alphabetDimension      =            theFilter->GetDimension(),
    siteCount            =            theFilter->GetPatternCount(),
    alphabetDimensionmod4  =          alphabetDimension-alphabetDimension%4;
    
    _CalcNode       *theNode               =            ((_CalcNode*) flatCLeaves (0));
    hyFloat  *   _hprestrict_ transitionMatrix
    =           theNode->GetCompExp(catID)->theData,
    result                 =            0.;
    
    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }
    
    for (long siteID = siteFrom; siteID < siteTo; siteID++) {
        hyFloat  *tMatrix = transitionMatrix,
        sum     = 0.;
        
        long siteState1 = lNodeFlags[siteOrdering.list_data[siteID]],
        siteState2 = lNodeFlags[siteCount + siteOrdering.list_data[siteID]];
        
        if (siteState1 >= 0)
            // a single character state; sweep down the appropriate column
        {
            if (siteState2 >= 0) { // both completely resolved;
                sum = tMatrix[siteState1*alphabetDimension + siteState2];
            } else { // first resolved, second is not
                hyFloat* childVector = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;
                tMatrix   +=  siteState1*alphabetDimension;
                if (alphabetDimension == 4) { // special case for nuc data
                    sum = tMatrix[0] * childVector[0] +
                    tMatrix[1] * childVector[1] +
                    tMatrix[2] * childVector[2] +
                    tMatrix[3] * childVector[3];
                } else {
                    for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
                        sum +=  tMatrix[c]   * childVector[c] +
                        tMatrix[c+1] * childVector[c+1] +
                        tMatrix[c+2] * childVector[c+2] +
                        tMatrix[c+3] * childVector[c+3];
                    
                    for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                        sum +=  tMatrix[c] * childVector[c];
                    }
                }
            }
            sum *= theProbs[siteState1];
        } else {
            if (siteState2 >=0 ) { // second resolved, but not the first
                hyFloat* childVector = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension;
                tMatrix                +=  siteState2;
                if (alphabetDimension == 4) { // special case for nuc data
                    sum = tMatrix[0] * childVector[0]  * theProbs[0]+
                    tMatrix[4] * childVector[1]  * theProbs[1]+
                    tMatrix[8] * childVector[2]  * theProbs[2]+
                    tMatrix[12] * childVector[3] * theProbs[3];
                    
                } else {
                    for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4*alphabetDimension) // 4 - unroll the loop
                        sum +=  tMatrix[0]   * childVector[c] * theProbs[c]+
                        tMatrix[alphabetDimension] * childVector[c+1] * theProbs[c+1]+
                        tMatrix[alphabetDimension+alphabetDimension] * childVector[c+2] * theProbs[c+2]+
                        tMatrix[alphabetDimension+alphabetDimension+alphabetDimension] * childVector[c+3] * theProbs[c+3];
                    
                    for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix += alphabetDimension) {
                        sum +=  tMatrix[0] * childVector[c] * theProbs[c];
                    }
                }
            } else
                // both unresolved
            {
                hyFloat *childVector1 = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension,
                *childVector2 = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;
                
                if (alphabetDimension == 4) { // special case for nuc data
                    sum = (tMatrix[0] * childVector2[0] + tMatrix[1] * childVector2[1] + tMatrix[2] * childVector2[2] + tMatrix[3] * childVector2[3])     * childVector1[0] * theProbs[0]+
                    (tMatrix[4] * childVector2[0] + tMatrix[5] * childVector2[1] + tMatrix[6] * childVector2[2] + tMatrix[7] * childVector2[3])     * childVector1[1] * theProbs[1]+
                    (tMatrix[8] * childVector2[0] + tMatrix[9] * childVector2[1] + tMatrix[10] * childVector2[2] + tMatrix[11] * childVector2[3])   * childVector1[2] * theProbs[2] +
                    (tMatrix[12] * childVector2[0] + tMatrix[13] * childVector2[1] + tMatrix[14] * childVector2[2] + tMatrix[15] * childVector2[3]) * childVector1[3] * theProbs[3];
                    
                } else {
                    for (long r = 0; r < alphabetDimension; r++) { // 4 - unroll the loop
                        if (childVector1[r] > 0.0) {
                            hyFloat sum2 = 0.0;
                            for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4) // 4 - unroll the loop
                                sum2   +=  tMatrix[0] * childVector2[c] +
                                tMatrix[1] * childVector2[c+1]+
                                tMatrix[2] * childVector2[c+2]+
                                tMatrix[3] * childVector2[c+3];
                            
                            for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix ++) {
                                sum2 +=  tMatrix[0] * childVector2[c];
                            }
                            
                            sum += sum2 * childVector1[r] * theProbs[r];
                        } else {
                            tMatrix += alphabetDimension;
                        }
                    }
                }
            }
        }
        if (storageVec) {
            storageVec [siteOrdering.list_data[siteID]] = sum;
        } else {
            if (sum <= 0.0) {
                return -INFINITY;
            } else {
                //printf ("%d: %g\n", siteID, sum);
                result += log(sum) * theFilter->theFrequencies.get (siteOrdering.list_data[siteID]);
            }
        }
    }
    
    return result;
}



//_______________________________________________________________________________________________

void     _TheTree::SampleAncestorsBySequence (_DataSetFilter const* dsf, _SimpleList const& siteOrdering, node<long>* currentNode, _AVLListX const* nodeToIndex, hyFloat const* iNodeCache,
                                              _List& result, _SimpleList* parentStates, _List& expandedSiteMap, hyFloat const* catAssignments, long catCount)

// must be called initially with the root node


// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual pattern ordering
// currentNode:                 the node index to sample for
// nodeToIndex:                 an AVL that maps the address of an internal node pointed to by node<long> to its order in the tree postorder traversal
// iNodeCache:                  internal node likelihood caches
// results:                     the list that will store sampled strings
// parentStates:                sampled states for the parent of the current node
// expandedSiteMap:             a list of simple lists giving site indices for each unique column pattern in the alignment
// catAssignments:              a vector assigning a (partition specific) rate category to each site (nil if no rate variation)
// catCount:                    the number of rate classes

// this needs to be updated to deal with traversal caches!
{
    long                      childrenCount     = currentNode->get_num_nodes();
    
    if (childrenCount) {
        long              siteCount                       = dsf->GetPatternCount  (),
        alphabetDimension              = dsf->GetDimension         (),
        nodeIndex                       = nodeToIndex->GetXtra (nodeToIndex->Find ((BaseRef)currentNode)),
        unitLength                     = dsf->GetUnitLength(),
        catBlockShifter                    = catAssignments?(dsf->GetPatternCount()*GetINodeCount()):0;
        
        
        _CalcNode *     currentTreeNode = ((_CalcNode*) flatTree (nodeIndex));
        _SimpleList     sampledStates     (dsf->GetSiteCountInUnits (), 0, 0);
        
        hyFloat  const *  transitionMatrix = (catAssignments|| !parentStates)?nil:currentTreeNode->GetCompExp()->theData;
        hyFloat  const *  conditionals     = catAssignments?nil:(iNodeCache + nodeIndex  * siteCount * alphabetDimension);
        hyFloat        *  cache            = new hyFloat [alphabetDimension];
        
        for (long           pattern = 0; pattern < siteCount; pattern++) {
            _SimpleList*    patternMap = (_SimpleList*) expandedSiteMap (siteOrdering.list_data[pattern]);
            if (catAssignments) {
                long localCatID = catAssignments[siteOrdering.list_data[pattern]];
                if (parentStates) {
                    transitionMatrix = currentTreeNode->GetCompExp(localCatID)->theData;
                }
                
                conditionals     = iNodeCache + localCatID*alphabetDimension*catBlockShifter + (pattern + nodeIndex  * siteCount) * alphabetDimension;
            }
            
            for (long site = 0; site < patternMap->lLength; site++) {
                long        siteID =   patternMap->list_data[site];
                
                hyFloat  randVal  = genrand_real2(),
                totalSum = 0.;
                
                hyFloat  const *   matrixRow;
                
                if  (parentStates == nil) {
                    matrixRow = theProbs;
                } else {
                    matrixRow = transitionMatrix + parentStates->list_data[siteID] * alphabetDimension;
                }
                
                for (long i = 0; i<alphabetDimension; i++) {
                    totalSum += (cache[i] = matrixRow[i]*conditionals[i]);
                }
                
                randVal *= totalSum;
                totalSum    = 0.0;
                long        sampledChar = -1;
                while       (totalSum < randVal) {
                    sampledChar ++;
                    totalSum += cache[sampledChar];
                }
                
                sampledStates.list_data[siteID] = sampledChar;
            }
            
            if (catAssignments == nil) {
                conditionals += alphabetDimension;
            }
        }
        
        delete [] cache;
        
        _SimpleList  conversion;
        _AVLListXL   conversionAVL (&conversion);
        
        _StringBuffer * sampledSequence = new _StringBuffer (siteCount*unitLength);
        _String  letterValue ((unsigned long) unitLength);
        for (long charIndexer = 0; charIndexer < sampledStates.countitems(); charIndexer++) {
            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(sampledStates.list_data[charIndexer]), unitLength, letterValue, &conversionAVL);
            (*sampledSequence) << letterValue;
        }
        sampledSequence->TrimSpace();
        result.AppendNewInstance(sampledSequence);
        //printf ("%d: %s\n", nodeIndex, sampledSequence->sData);
        
        for (long child = 1L; child <= childrenCount; child ++) {
            SampleAncestorsBySequence (dsf, siteOrdering, currentNode->go_down(child), nodeToIndex, iNodeCache, result, &sampledStates, expandedSiteMap, catAssignments, catCount);
        }
    }
}

//_______________________________________________________________________________________________

_List*   _TheTree::RecoverAncestralSequences (_DataSetFilter const* dsf,
                                              _SimpleList const& siteOrdering,
                                              _List const& expandedSiteMap,
                                              hyFloat * iNodeCache,
                                              hyFloat const* catAssignments,
                                              long catCount,
                                              long* lNodeFlags,
                                              _Vector * lNodeResolutions,
                                              bool              alsoDoLeaves
                                              )


// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual pattern ordering
// expandedSiteMap:             a list of simple lists giving site indices for each unique column pattern in the alignment
// iNodeCache:                  internal node likelihood caches
// catAssignments:              a vector assigning a (partition specific) rate category to each site
// catCount:                    the number of rate classes
// alsoDoLeaves:                if true, also return ML reconstruction of observed (or partially observed) sequences
{
    long            patternCount                     = dsf->GetPatternCount  (),
    alphabetDimension                = dsf->GetDimension         (),
    unitLength                       = dsf->GetUnitLength        (),
    iNodeCount                       = GetINodeCount             (),
    leafCount                        = GetLeafCount              (),
    siteCount                        = dsf->GetSiteCountInUnits    (),
    allNodeCount                     = 0,
    stateCacheDim                    = (alsoDoLeaves? (iNodeCount + leafCount): (iNodeCount));
    
    long            *stateCache                     = new long [patternCount*(iNodeCount-1)*alphabetDimension],
                    *leafBuffer                     = new long [(alsoDoLeaves?leafCount*patternCount:1)*alphabetDimension],
                    *initiaStateCache               = stateCache;
    
    // a Patterns x Int-Nodes x CharStates integer table
    // with the best character assignment for node i given that its parent state is j for a given site
    
    hyFloat          *buffer                         = new hyFloat [alphabetDimension];
    // iNodeCache will be OVERWRITTEN with conditional pair (i,j) conditional likelihoods
    
    
    _SimpleList     taggedInternals (iNodeCount, 0, 0),
    postToIn;
    
    MapPostOrderToInOrderTraversal (postToIn);
    // all nodes except the root
    
    allNodeCount = iNodeCount + leafCount - 1;
    
    if (iNodeCache) {
        for  (long nodeID = 0; nodeID < allNodeCount; nodeID++) {
            long    parent_index = flatParents.get (nodeID),
                    node_index   = nodeID;
            
            bool    is_leaf     = nodeID < flatLeaves.countitems();
            
            
            if (!is_leaf) {
                node_index -=  flatLeaves.countitems();
                AddBranchToForcedRecomputeList (node_index);
            }
            
            hyFloat * parentConditionals = iNodeCache + parent_index * alphabetDimension * patternCount;
            
            if (taggedInternals.get(parent_index) == 0L) {
                // mark the parent for update and clear its conditionals if needed
                taggedInternals[parent_index]     = 1L;
                InitializeArray(parentConditionals, patternCount*alphabetDimension, 1.);
            }
            
            _CalcNode *          tree_node_object = is_leaf? ((_CalcNode*) flatCLeaves (node_index)):((_CalcNode*) flatTree    (node_index));
            hyFloat  const*      transition_matrix = nil;
            
            if (!catAssignments) {
                _Matrix* comp_exp = tree_node_object->GetCompExp();
                if (!comp_exp) {
                    hy_global::HandleApplicationError(_String ("Internal error in ") & __PRETTY_FUNCTION__ & ". Transition matrix not computed for " & *tree_node_object->GetName());
                    delete [] stateCache; delete [] leafBuffer; delete [] buffer;
                    return nil;
                }
                transition_matrix = comp_exp->theData;
            }
            
            // this will need to be toggled on a per site basis
            hyFloat  *       childVector;
            
            if (!is_leaf) {
                childVector = iNodeCache + (node_index * patternCount) * alphabetDimension;
            }
            
            for (long siteID = 0; siteID < patternCount; siteID++, parentConditionals += alphabetDimension) {
                if (catAssignments) {
                    transition_matrix = tree_node_object->GetCompExp(catAssignments[siteOrdering.list_data[siteID]])->theData;
                }
                
                hyFloat  const *tMatrix = transition_matrix;
                if (is_leaf) {
                    long siteState = lNodeFlags[node_index*patternCount + siteOrdering.list_data[siteID]] ;
                    if (siteState >= 0L) { // a fully resolved leaf
                        tMatrix  +=  siteState;
                        for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                            parentConditionals[k] *= *tMatrix;
                        }
                        if (alsoDoLeaves) {
                            InitializeArray (leafBuffer, alphabetDimension, (const long)siteState);
                            leafBuffer += alphabetDimension;
                        }
                        
                        continue;
                    } else {// an ambiguous leaf
                        childVector = lNodeResolutions->theData + (-siteState-1L) * alphabetDimension;
                    }
                    
                }
                
                // now repopulate this vector as necessary -- if we are here this means
                // that the subtree below has been completely processed,
                // the i-th cell of childVector contains the likelihood of the _optimal_
                // assignment in the subtree below given that the character at the current
                // node is i.
                
                // hence, given parent state 'p', we optimize
                // max_i pr (p->i) childVector [i] and store it in the p cell of vector childVector
                
                hyFloat overallMax                     = 0.0;
                
                long       *stateBuffer                   = is_leaf?leafBuffer:stateCache;
                
                // check for degeneracy
                
                bool completely_unresolved = ArrayAll (childVector, alphabetDimension, [] (hyFloat x, unsigned long) {return x == 1.;});
                
                if (completely_unresolved) {
                    InitializeArray(stateBuffer, alphabetDimension, -1L);
                } else {
                    for (long p = 0L; p < alphabetDimension; p++) {
                        hyFloat max_lik = 0.;
                        long       max_idx = 0L;
                        
                        for (long c = 0L; c < alphabetDimension; c++) {
                            hyFloat thisV = tMatrix[c] * childVector[c];
                            if (thisV > max_lik) {
                                max_lik = thisV;
                                max_idx = c;
                            }
                        }
                        
                        stateBuffer [p] = max_idx;
                        buffer [p]      = max_lik;
                        
                        if (max_lik > overallMax) {
                            overallMax = max_lik;
                        }
                        
                        tMatrix += alphabetDimension;
                    }
                    
                    if (overallMax > 0.0 && overallMax < _lfScalingFactorThreshold) {
                        for (long k = 0L; k < alphabetDimension; k++) {
                            buffer[k] *= _lfScalerUpwards;
                        }
                    }
                    
                    // buffer[p] now contains the maximum likelihood of the tree
                    // from this point forward given that parent state is p
                    // and stateBuffer[p] stores the maximizing assignment
                    // for this node
                    
                    for (long k = 0; k < alphabetDimension; k++) {
                        if (stateBuffer[k] >= 0L) {
                            parentConditionals[k] *= buffer[k];
                        }
                    }
                }
                
                if (is_leaf) {
                    if (alsoDoLeaves) {
                        leafBuffer += alphabetDimension;
                    }
                } else {
                    stateCache += alphabetDimension;
                }
                
                childVector += alphabetDimension;
            }
        }
    }
    
    _List      *result = new _List;
    
    if (!iNodeCache) { // special case for 2 sequence data; root state is actually the 1st of the two sequences
        //result->AppendNewInstance
        result->AppendNewInstance (dsf->GetSequenceCharacters(0L));
    } else {
        for (long k = 0; k < stateCacheDim; k++) {
            result->AppendNewInstance (new _String((unsigned long)siteCount*unitLength));
        }
     
        hyFloat * _hprestrict_ rootConditionals = iNodeCache + alphabetDimension * ((iNodeCount-1)  * patternCount);
        _SimpleList  parentStates (stateCacheDim,0,0),
        conversion;
        
        stateCache -= patternCount*(iNodeCount-1)*alphabetDimension;
        if (alsoDoLeaves) {
            leafBuffer -= patternCount*leafCount*alphabetDimension;
        }
        
        _AVLListXL    conversionAVL (&conversion);
        _String       codeBuffer    ((unsigned long)unitLength);
        
        for (long siteID = 0; siteID < patternCount; siteID++, rootConditionals += alphabetDimension) {
            hyFloat max_lik = 0.;
            long       max_idx = 0;
            
            long howManyOnes = 0;
            for (long k = 0; k < alphabetDimension; k++) {
                howManyOnes += rootConditionals[k]==1.;
            }
            
            _SimpleList const*    patternMap = (_SimpleList const*) expandedSiteMap.GetItem(siteOrdering.list_data[siteID]);
            
            if (howManyOnes != alphabetDimension) {
                for (long c = 0; c < alphabetDimension; c++) {
                    hyFloat thisV = theProbs[c] * rootConditionals[c];
                    if (thisV > max_lik) {
                        max_lik = thisV;
                        max_idx = c;
                    }
                }
                
                parentStates.list_data[iNodeCount-1] = max_idx;
                for  (long nodeID = iNodeCount-2; nodeID >=0 ; nodeID--) {
                    long parentState = parentStates.list_data[flatParents.list_data [nodeID+flatLeaves.lLength]];
                    if (parentState == -1) {
                        parentStates.list_data[nodeID] = -1;
                    } else {
                        parentStates.list_data[nodeID] = stateCache[(patternCount*nodeID+siteID)*alphabetDimension + parentState];
                    }
                }
                if (alsoDoLeaves)
                    for  (long nodeID = 0; nodeID <leafCount ; nodeID++) {
                        long parentState = parentStates.list_data[flatParents.list_data [nodeID]];
                        if (parentState == -1) {
                            parentStates.list_data[nodeID+iNodeCount] = -1;
                        } else {
                            parentStates.list_data[nodeID+iNodeCount] = leafBuffer[(patternCount*nodeID+siteID)*alphabetDimension + parentState];
                        }
                    }
            } else {
                parentStates.Populate(stateCacheDim,-1,0);
            }
            
            for  (long nodeID = 0; nodeID < stateCacheDim ; nodeID++) {
                dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(parentStates.list_data[nodeID]), unitLength, codeBuffer, &conversionAVL);
                _String  *sequence   = (_String*) (*result)(nodeID<iNodeCount?postToIn.list_data[nodeID]:nodeID);
                
                for (long site = 0; site < patternMap->lLength; site++) {
                    unsigned long offset = patternMap->list_data[site]*unitLength;
                    for (long charS = 0; charS < unitLength; charS ++) {
                        sequence->set_char (offset + charS, codeBuffer.char_at(charS));
                    }
                }
            }
        }
    }
    
    delete [] initiaStateCache;
    delete [] leafBuffer;
    delete [] buffer;
    
    return result;
}

//_______________________________________________________________________________________________

void     _TheTree::SetupCategoryMapsForNodes (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        iterator->SetupCategoryMap (containerVariables,classCounter,multipliers);
    }
}



//_______________________________________________________________________________________________

hyFloat   _TheTree::Process3TaxonNumericFilter (_DataSetFilterNumeric* dsf, long catID)
{
    
    hyFloat *l0 =  dsf->probabilityVectors.theData +
    dsf->categoryShifter * catID + dsf->theNodeMap.list_data[0]*dsf->shifter,
    *l1 = dsf->probabilityVectors.theData +
    dsf->categoryShifter * catID + dsf->theNodeMap.list_data[1]*dsf->shifter,
    *l2 = dsf->probabilityVectors.theData +
    dsf->categoryShifter * catID + dsf->theNodeMap.list_data[2]*dsf->shifter,
    * matrix0 = ((_CalcNode*)(LocateVar(theRoot->get_node(1)->in_object)))->GetCompExp(catID)->theData,
    * matrix1 = ((_CalcNode*)(LocateVar(theRoot->get_node(2)->in_object)))->GetCompExp(catID)->theData,
    * matrix2 = ((_CalcNode*)(LocateVar(theRoot->get_node(3)->in_object)))->GetCompExp(catID)->theData,
    overallResult = 0.;
    
    long        patternCount =  dsf->GetPatternCount();
    
    hyFloat  currentAccumulator = 1.;
    
    for (long patternIndex = 0; patternIndex < patternCount; patternIndex ++, l0+=4, l1+=4, l2+=4) {
        hyFloat rp0 = l0[0] * matrix0[0]+ l0[1]  * matrix0[1]  + l0[2] * matrix0[2]  + l0[3] * matrix0[3];
        hyFloat rp1 = l0[0] * matrix0[4]+ l0[1]  * matrix0[5]  + l0[2] * matrix0[6]  + l0[3] * matrix0[7];
        hyFloat rp2 = l0[0] * matrix0[8]+ l0[1]  * matrix0[9]  + l0[2] * matrix0[10] + l0[3] * matrix0[11];
        hyFloat rp3 = l0[0] * matrix0[12]+ l0[1] * matrix0[13] + l0[2] * matrix0[14] + l0[3] * matrix0[15];
        
        rp0 *= l1[0] * matrix1[0] + l1[1] * matrix1[1]  + l1[2] * matrix1[2]  + l1[3] * matrix1[3];
        rp1 *= l1[0] * matrix1[4] + l1[1] * matrix1[5]  + l1[2] * matrix1[6]  + l1[3] * matrix1[7];
        rp2 *= l1[0] * matrix1[8] + l1[1] * matrix1[9]  + l1[2] * matrix1[10] + l1[3] * matrix1[11];
        rp3 *= l1[0] * matrix1[12]+ l1[1] * matrix1[13] + l1[2] * matrix1[14] + l1[3] * matrix1[15];
        
        rp0 *= l2[0] * matrix2[0] + l2[1] * matrix2[1]  + l2[2] * matrix2[2]  + l2[3] * matrix2[3];
        rp1 *= l2[0] * matrix2[4] + l2[1] * matrix2[5]  + l2[2] * matrix2[6]  + l2[3] * matrix2[7];
        rp2 *= l2[0] * matrix2[8] + l2[1] * matrix2[9]  + l2[2] * matrix2[10] + l2[3] * matrix2[11];
        rp3 *= l2[0] * matrix2[12]+ l2[1] * matrix2[13] + l2[2] * matrix2[14] + l2[3] * matrix2[15];
        
        hyFloat  result = theProbs[0]*rp0+
        theProbs[1]*rp1+
        theProbs[2]*rp2+
        theProbs[3]*rp3;
        
        
        if (result<=0.0) {
            return -INFINITY;
        }
        
        long patternFreq = dsf->theFrequencies[patternIndex];
        for  (long freqIterator = 0; freqIterator < patternFreq; freqIterator++) {
            hyFloat tryMultiplication = currentAccumulator*result;
            if (tryMultiplication > 1.e-300) {
                currentAccumulator = tryMultiplication;
            } else {
                overallResult += myLog (currentAccumulator);
                currentAccumulator = result;
            }
        }
    }
    return overallResult + myLog (currentAccumulator);
}

//_______________________________________________________________________________________________

void    _TheTree::_RemoveNodeList (_SimpleList const& clean_indices) {
    if (compExp) {
        DeleteAndZeroObject(compExp);
    }
    clean_indices.Each([] (long var_idx, unsigned long) -> void {
        DeleteVariable(var_idx, true);
    });
}
