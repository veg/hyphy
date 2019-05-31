
ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "BayesianGraphicalModel";
}		


function runTest () {
   ASSERTION_BEHAVIOR = TRUE; /* print warning to console and go to the end of the execution list */
   testResult = FALSE;
  
  // TODO: further testing of the BayesianGraphicalModel function... 
  // currently just confirming that we can create a simple model.
  LoadFunctionLibrary ("libv3/tasks/bayesgraph.ibf");
  nodes = {};
  num_nodes = 10;
  max_parents = 2;
	for (k = 0; k < num_nodes; k += 1) {
	    node_name = ""+ (k + 1);
	    nodes + bgm.add_discrete_node (node_name, max_parents, 0, 2);
	}
  BayesianGraphicalModel gen_bgm = (nodes);

  tempFilePathBGM = './../../data/tempFileTesting-bgm-' + Random(0,1);
  fprintf (tempFilePathBGM, gen_bgm);
  fscanf (tempFilePathBGM, String, bgmPrintOut);
  assert(bgmPrintOut == "Log Likelihood =               0;", "Failed to initialize a bgm with a log likelihood of 0");
  
  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};
  Topology Topo1 = (a,b,c);
  Tree  Tree1 = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);

  assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = ('string');", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with a string");
  assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = (100);", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with a number");
  assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = (matrix1);", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with a matrix");
  assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = (Topo1);", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with a matrix");
  assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = (Tree1);", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with a matrix");
  
  // TODO: segfaults with a list
  //assert (runCommandWithSoftErrors ("BayesianGraphicalModel gen_bgm = (list1);", "n call to BGM = ... must evaluate to associative array"), "Failed error checking for calling bgm with the wrong type of associative list");
  
  
  testResult = 1;

  return testResult;
}
