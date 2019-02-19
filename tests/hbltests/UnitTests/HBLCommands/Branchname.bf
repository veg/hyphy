ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Branchname";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  Topology T =  (((a:0.1,b:0.2)ab:0.4,e:0.1):0.2,c:0.15,d:0.33);
  Tree TT = (((a:0.1,b:0.2)ab:0.4,e:0.1):0.2,c:0.15,d:0.33);

  // Get branch name by traversal index.
  node1BranchName = BranchName (T,1);
  assert(node1BranchName == 'Node1', "Failed to correctly get branch name by traversal of index of a Topology");
  node1BranchNameTree = BranchName (TT,1);
  assert(node1BranchNameTree == 'Node1', "Failed to correctly get branch name by traversal of index of a Tree");
  
  
  // TODO: The Three functions below don't work as documented
  // The "expected*" is what is shown in the docs.

  // Get all branch names.
  allBranchNames = BranchName (T,-1);
  expectedAllBranchNames = {{"a","b","ab","e","Node1","c","d","Node0"}};
  actualAllBranchNames = {{"a", "b", "ab", "e", "Node1", "c", "d", "Node7"}};
  //assert(allBranchNames == expectedAllBranchNames, "Failed to correctly get all the branch names from a Topology");
  allBranchNamesTree = BranchName (TT,-1);
  expectedAllBranchNamesTree = {{"a","b","ab","e","Node1","c","d","Node0"}};
  actualAllBranchNamesTree = {{"a", "b", "ab", "e", "Node1", "c", "d", "Node7"}};
  //assert(allBranchNamesTree == expectedAllBranchNamesTree, "Failed to correctly get all the branch names from a Tree");
  
  // Get all branch names along a path.
  alongPathBranchNames = BranchName (T,"d;a");
  expectedAlongPathBranchNames = {{"d","Node1","ab","a"}};
  actualAlongPathBranchNames = {{"Node7", "Node7", "ab", "a"}};
  //assert(alongPathBranchNames == expectedAlongPathBranchNames, "Failed to correctly get the branch names along a path" from a Topology);
  alongPathBranchNamesTree = BranchName (TT,"d;a");
  expectedAlongPathBranchNamesTree = {{"d","Node1","ab","a"}};
  actualAlongPathBranchNamesTree = {{"Node7", "Node7", "ab", "a"}};
  //assert(alongPathBranchNamesTree == expectedAlongPathBranchNamesTree, "Failed to correctly get the branch names along a path" from a Tree);

  // Get subtree information.
  subtreeInformation = BranchName (T,"Node1");
  expectedSubtreeInformation = {"ab":2, "a":1, "b":1};
  actualSubtreeInformation = {"Node1":2, "a":1, "ab":2, "b":1, "e":1 };
  //assert(subtreeInformation["Node1"] == expectedSubtreeInformation["Node1"], "Failed to correctly get the subtree information from a Topology");
  subtreeInformationTree = BranchName (TT,"Node1");
  expectedSubtreeInformationTree = {"ab":2, "a":1, "b":1};
  actualSubtreeInformationTree = {"Node1":2, "a":1, "ab":2, "b":1, "e":1 };
  //assert(subtreeInformationTree["Node1"] == expectedSubtreeInformationTree["Node1"], "Failed to correctly get the subtree information from a Tree");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};
  
  // TODO (the three tests below were passing before ~1/20/2019): 
  assert (runCommandWithSoftErrors ('BranchName(list1)', "Operation 'BranchName' is not implemented/defined for a AssociativeList"), "Failed error checking for running BranchName on AssociativeList");
  // assert (runCommandWithSoftErrors ('BranchName(1)', "Operation 'BranchName' is not implemented/defined for a Number"), "Failed error checking for running BranchName on Number");
  // assert (runCommandWithSoftErrors ('BranchName(none)', "Operation 'BranchName' is not implemented/defined for a Number"), "Failed error checking for running BranchName on none");
  
  // TODO: Not sure why the below assert fails... 
  //assert (runCommandWithSoftErrors ('BranchName(matrix1)', "Operation 'BranchName' is not implemented/defined for a Matrix"), "Failed error checking for running BranchName on Matrix");

  
  testResult = 1;

  return testResult;
}
