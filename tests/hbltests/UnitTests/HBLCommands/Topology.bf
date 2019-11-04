ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Topology";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
 
    ACCEPT_ROOTED_TREES = 0;
   Topology single_leaf = (a,b);
   fprintf (stdout, Format (single_leaf,0,0), "\n"); 
   

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Many of the functions that can be performed on Topologys (branchcount, branchlength, etc.) have their own unit tests.
  // This unit test is focused on the Topology data structure apptempting to cover things not addressed by other tests.
  
  // Rooted vs. unrooted is handeled by the `ACCEPT_ROOTED_TopologyS` environmental variable
  Topology simpleTopology = ((a,b),(c,d));
  ACCEPT_ROOTED_TREES = 0;
  Topology simpleTopologyUnrooted = ((a,b),(c,d));
  ACCEPT_ROOTED_TREES = 1;
  Topology simpleTopologyRooted = ((a,b),(c,d));
  assert(BranchCount(simpleTopology) == BranchCount(simpleTopologyUnrooted), "Failed to create an unrooted Topology by default, i.e. if ACCEPT_ROOTED_TopologyS is not defined");
  assert((BranchCount(simpleTopologyRooted) - BranchCount(simpleTopologyUnrooted)) == 1, "A rooted version of a Topology did not have one more internal branch than the unrooted version");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};
  assert (runCommandWithSoftErrors ('Topology listTopology = list1', "Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create Topology with list"); 
  assert (runCommandWithSoftErrors ('Topology matrixTopology = matrix1', "Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create Topology with matrix"); 
  assert (runCommandWithSoftErrors ('Topology numberTopology = 3', "Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology"), "Failed error checking for trying to create Topology with number"); 

  ACCEPT_ROOTED_TREES = 1;
  Topology twoNodes = (a),b;
  ACCEPT_ROOTED_TREES = 0;

  // TODO: Tree and Topology behave differently... Tree throws a "cannot construct emtpy trees" error wheras Topology does not... I don't think this is desired  
  //Topology twoNodes = (a),b;
  //assert (runCommandWithSoftErrors ('Topology twoNodes = (a),b;', "Cannot constuct empty Topologys"), "Failed error checking for trying to construct a Topology with only two nodes");
  
    
  testResult = 1;

  return testResult;
}
