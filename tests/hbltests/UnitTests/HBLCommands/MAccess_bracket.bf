ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
return runATest ();


function getTestName () {
  return "MAccess_bracket";
}		

function dict_callback_action (key, value) {
    sum += Abs (value);
}

function dict_callback_filter (key) {
    return key != "key3";
}

function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Element access operation for matrices, trees and strings.

  // String.
  exampleString = "abcd";

  exampleStringSecond = exampleString[1];
  exampleStringSecondOperationInBracket = exampleString[3-1+1];
  exampleStringOutOfRange = exampleString[4];
  exampleStringNegative = exampleString[-1];

  assert(exampleStringSecond == "b", "Failed element access on string");
  assert(exampleStringSecondOperationInBracket == "d", "Failed element access on string with mathematical operation in bracket");
  assert(exampleStringOutOfRange == "", "Failed to return empty string for out of range element access on string (4)");
  assert(exampleStringNegative == "d", "Failed to last element for negative access (-1)");

  // Matrix. Accesses as expected, if only one index is given to a 2d array, the array is unwrapped prior to access.
  exampleMatrix = {{1,2,3}{4,5,6}{7,8,9}};
  
  exampleMatrixFirstThird = exampleMatrix[0][2];
  exampleMatrixOnlySecond = exampleMatrix[1];
  exampleMatrixOnlyFourth = exampleMatrix[3];
  
  assert(exampleMatrixFirstThird == 3, "Failed simple element access on 2d array");
  assert(exampleMatrixOnlySecond == 2, "Failed one argument element access on 2d array (i.e. M=3x3Matrix M[2])");
  assert(exampleMatrixOnlyFourth == 4, "Failed one argument element access on 2d array (i.e. M=3x3Matrix M[3])");

  // Tree. Tree[i] returns a string representing the subtree eminating from the ith node.
  Tree exampleTree = ((1,2),(3,4),5);

  exampleTreeFirst = exampleTree[0];
  exampleTreeSecond = exampleTree[1];

  assert(exampleTreeFirst == "(1,2)Node1", "Failed element access on a tree");
  assert(exampleTreeSecond =="(3,4)Node4", "Failed element acess on a tree");

  // Topology (same as tree). Topology[i] returns a string representing the subtopology eminating from the ith node.
  Topology exampleTopology = ((1,2),(3,4),5);

  exampleTopologyFirst = exampleTopology[0];
  exampleTopologySecond = exampleTopology[1];

  assert(exampleTopologyFirst == "(1,2)Node1", "Failed element access on a topology");
  assert(exampleTopologySecond == "(3,4)Node4", "Failed element acess on a topology");

  // Associative list. Access via the quoted key.
  exampleList = {"key1": "val1", "key2": "val2", "key3": {"subkey1": "subvalue1", "subkey2": "subvalue2"}};
  assert(exampleList["key1"] == "val1", "Failed to access a value from an associative list with the `list['key']` syntax");
  // TODO: Should nested associative lists be accessable without having to wrap the first call in parenthesis? See below.
  assert((exampleList["key3"])["subkey1"] == "subvalue1", "Failed to access a value from a nested associative list when using parentesis `(list['key1'])[key2']` synatax");
  
  
  //---------------------------------------------------------------------------------------------------------
  // callback on associative lists
  //---------------------------------------------------------------------------------------------------------

  sum = 0;
  exampleList["dict_callback_action"][""];
  assert(sum == 10, 'Incorrect result using dict callback exampleList["dict_callback_action"]');
   
  sum = 0;
  exampleList["dict_callback_action"]["dict_callback_filter"];
  assert(sum == 8, 'Incorrect result using dict callback exampleList["dict_callback_action"]["dict_callback_filter"]');
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  
  // TODO: error handling for trying element access lists or numbers should be improved.
  
  assert (runCommandWithSoftErrors ('exampleList["getTestName"][""]', 'The first argument in an iterator call for Associative Arrays must be a valid identifier of a function taking two arguments \(key, value\)'), "Failed to return the expected error when providing an invalid callback");
  
  assert (runCommandWithSoftErrors ('exampleList["dict_callback_action"]["getTestName"]', 'The second argument in an iterator call for Associative Arrays must be either empty or a valid identifier of a function taking a single argument'), "Failed to return the expected error when providing an invalid callback");

  // Test for current behaviour for list (return 0 when trying element access on list)
  assert(exampleList[0] == 0, "Failed to return zero when trying element access on list");
  assert(exampleList[100] == 0, "Failed to return zero when trying out of range element access on lilst");


 

  testResult = 1;

  return testResult;
}

