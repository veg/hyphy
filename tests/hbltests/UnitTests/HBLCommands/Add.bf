ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "+";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = FALSE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Add numbers 
  assert(1+2 == 3, "Failed to add two integers");
  assert(1.00001+3.33333 == 4.33334, "Failed to add two floats");
  assert(-100+1 == -99, "Failed to add a positive int to a negative int");
  // Add Matrices
  assert({{1,2,3}} + {{3,2,1}} == {{4,4,4}}, "Failed to add two 1d matrices");
  assert({{1,2}{3,4}} + {{1,2}{3,4}} == {{2,4}{6,8}}, "Failed to add two 2d matrices");
  // Concatenate strings
  assert("juxta"+"position" == "juxtaposition", "Failed to concatenate two strings");
  // Add two associative arrays (second elemet as a new key-value pair to the first element with the key of "1")
  list1 = {"key": "value"};
  list2 = {"key2": "value2", "key3": "value3"};
  expectedCombinedList = {"key": "value", "1": {"key2": "value2", "key3": "value3"}};
  list1 + list2;
  list1Subset = list1["1"];
  list1SubsetVal2 = list1Subset["key2"];
  assert(list1SubsetVal2 == "value2", "Failed two add two associative arrays");
  // Adding different types (string and numbers)
  assert("test" + 1 == "test1", "Failed to add a string and number (in that order); should convert the number to a string and concatenate it to the end of the string");
  assert("test" + 1 + 2 == "test12", "Failed to add a string and two numbers (in that order); should just return the sum of the two numbers");
  assert(1 + "test" == 1, "Failed to add a number and a string (in that order); should just return the number");
  assert(1 + "4" + 2 == 7, "Failed to add a number, a string (which is converted to a number) and another number (in that order; i.e. number first); should return the sum of the numbers");


  // test recursive sum on objects
  nested = {
    "0" : 1,
    "1" : {{2,3,4,5}},
    "2" : {
        "A" : {{6,7,8}},
        "B" : {"0" : 9, "1" : 10}
    }
  };
  
  
  sum = 0;
  for (i = 0; i < 10; i+=1) {
    sum += (+nested);
  }
  
  assert(sum == 550, "Failed to correctly sum up a nested object");
  

  Topology T = ((a,b)N1,c,d);
  // + on topologies works in place
  T + {"NAME":"e", "WHERE": "b", "PARENT": "f"};
  Topology expectedNewTopology = ((a,(b,e)f)N1,c,d);
  
    
  assert((T == expectedNewTopology) == " ", "Failed to add a node to a topology (1)");
  T + {"NAME":"Z", "WHERE": "f", "PARENT": "N2"};
  Topology expectedNewTopology = ((a,((b,e)f,Z)N2)N1,c,d);
  assert((T == expectedNewTopology) == " ", "Failed to add a node to a topology (2)");

  T + {"NAME":"Z", "WHERE": "f"};
  Topology expectedNewTopology = ((a,((b,e,Z)f,Z)N2)N1,c,d);
  assert((T == expectedNewTopology) == " ", "Failed to add a node to a topology without creating a new parent (internal node): " + T + " != " + expectedNewTopology);

  T + {"NAME":"X", "WHERE": "Z"};
  Topology expectedNewTopology = ((a,((b,e,(X)Z)f,Z)N2)N1,c,d);
  assert((T == expectedNewTopology) == " ", "Failed to add a node to a topology without creating a new parent (terminal node): " + T + " != " + expectedNewTopology);

  T + {"PARENT":"XZ", "WHERE": "Z"};
  Topology expectedNewTopology = ((a,((b,e,((X)Z)XZ)f,Z)N2)N1,c,d);
  assert((T == expectedNewTopology) == " ", "Failed to add a node to a topology without creating a new parent (terminal node): " + T + " != " + expectedNewTopology);

  assert (runCommandWithSoftErrors ('T + {"NAME":"Z", "WHERE": "missing", "PARENT": "N2"}', "Attachment node must be an exiting non-root node "), "Failed error adding a node to a non-existing node");
  assert (runCommandWithSoftErrors ('T + {"NAME":"Z"}', "Missing/invalid mandatory argument"), "Failed error checking when not specifying all required arguments");
  assert (runCommandWithSoftErrors ('T + {"WHERE":"X"}', "Either 'NAME' or 'PARENT'"), "Failed error checking when not specifying all required arguments");
 
  Tree T = ((a:0.1,b:0.2)N1:0.3,c:0.1,d:0.2);
  T + {"NAME":"e", "WHERE": "b", "LENGTH" : 0.05, "PARENT_LENGTH" : 0.01, "PARENT": "p_of_b"};
  assert(BranchLength (T, "e") == 0.05, "Incorrectly set branch length for an added branch");
  assert(BranchLength (T, "p_of_b") == 0.01, "Incorrectly set branch length for the parent of an added branch");
  T + {"NAME":"k", "WHERE": "b", "PARENT_LENGTH" : 0.011, "PARENT": "p_of_k"};
  assert(BranchLength (T, "k") == -1, "Incorrectly set missing branch length value for added branch");
  assert(BranchLength (T, "p_of_k") == 0.011, "Incorrectly set branch length for the parent of an added branch [no new branch length specified]");
  T + {"NAME":"n", "WHERE": "k", "LENGTH" : 0.012, "PARENT": "p_of_n"};
  assert(BranchLength (T, "p_of_n") == -1, "Incorrectly set missing branch length value for the parent of an added branch");
  assert(BranchLength (T, "n") == 0.012, "Incorrectly set branch length for an added branch [no parent branch length specified]");
  


  // TODO: Adding different types (string and matrix)
  stringPlusMatrix = "string+" + {{1,2}};
  assert(stringPlusMatrix == 
"string+{
{1, 2} 
}", "Failed to add a string and a matrix (in that order); should convert the matrix to a string and concatenate it to the end of the string");
   
 
  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  oneDMatrix = {{1,2}};
  twoDMatrix = {{3,4}{5,6}};

  assert (runCommandWithSoftErrors ('T2+T2', "An invalid argument"), "Failed error checking for trying to add two topologies (you can only add a node, i.e. associative array, to a topology)");
  assert (runCommandWithSoftErrors ('TT+TT', "An invalid argument"), "Failed error checking for trying to add trees");
  assert (runCommandWithSoftErrors ('oneDMatrix + twoDMatrix', "Incompatible dimensions when trying to add or subtract matrices"), "Failed error checking for adding matrices of different dimensions");
  

  testResult = TRUE;

  return testResult;
}
