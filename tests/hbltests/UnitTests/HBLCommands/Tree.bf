ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Tree";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Many of the functions that can be performed on trees (branchcount, branchlength, etc.) have their own unit tests.
  // This unit test is focused on the tree data structure apptempting to cover things not addressed by other tests.
  
  // Rooted vs. unrooted is handeled by the `ACCEPT_ROOTED_TREES` environmental variable
  Tree simpleTree = ((a,b),(c,d));
  ACCEPT_ROOTED_TREES = 0;
  Tree simpleTreeUnrooted = ((a,b),(c,d));
  ACCEPT_ROOTED_TREES = 1;
  Tree simpleTreeRooted = ((a,b),(c,d));
  assert(BranchCount(simpleTree) == BranchCount(simpleTreeUnrooted), "Failed to create an unrooted tree by default, i.e. if ACCEPT_ROOTED_TREES is not defined");
  assert((BranchCount(simpleTreeRooted) - BranchCount(simpleTreeUnrooted)) == 1, "A rooted version of a tree did not have one more internal branch than the unrooted version");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};
  assert (runCommandWithSoftErrors ('Tree listTree = list1', "Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create tree with list"); 
  assert (runCommandWithSoftErrors ('Tree matrixTree = matrix1', "Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create tree with matrix"); 
  assert (runCommandWithSoftErrors ('Tree numberTree = 3', "Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create tree with number"); 


  ACCEPT_ROOTED_TREES = 1;
  
  Tree twoNodes = (a),b;
  ACCEPT_ROOTED_TREES = 0;
  assert (runCommandWithSoftErrors ('Tree twoNodes = (a),b;', "Cannot consrtuct empty trees"), "Failed error checking for trying to construct a tree with only two nodes");
   
  testResult = 1;

  return testResult;
}
