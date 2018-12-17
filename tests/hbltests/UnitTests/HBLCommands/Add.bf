ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "+";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

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
  assert(1 + "test" + 2 == 3, "Failed to add a number, a string and another number (in that order; i.e. number first); should return the sum of the numbers");


  // Can't get the two tests below to work
  // Add topologies
  //Topology T = ((a,b)N1,c,d);
  //t = T + {"NAME":"e", "WHERE": "b", "PARENT": "f"};
  //fprintf (stdout, "t: ", t, "\n");
  //fprintf (stdout, "T: ", T, "\n");
  //expectedNewTopology = ((a,(b,e)f)N1,c,d);
  //assert(T == expectedNewTopology, "Failed to add a node to a topology");

  // Adding different types (string and matrix)
  //stringPlusMatrix = "string" + {{1,2}};
  //fprintf (stdout, "stringPlusMatrix: ", stringPlusMatrix, "\n");
  //assert(stringPlusMatrix == "test{{1,2}}", "Failed to add a string and a matrix (in that order); should convert the matrix to a string and concatenate it to the end of the string");
   
 
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
  

  testResult = 1;

  return testResult;
}
