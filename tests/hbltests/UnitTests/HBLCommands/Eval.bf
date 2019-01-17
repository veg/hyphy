ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Eval";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Evaluate a function in a string.
  assert(Eval("3+3*13") == 42, "Failed to evaluate a standard function");
  assert(Eval("3^3+2") == 29, "Failed to evaluate a function with exponents");
  assert(Eval("(3+3)*13") == 78, "Failed to evaluate a function with parentheses");
  assert(Eval("1+1+string") == 2, "Failed to evaluate a function with string as one of the arguments");
  
  // Other accepted types
  assert(Eval(none) == 0, "Failed to return zero when evaluation none");
  assert(Eval(1+1) == 2, "Failed to evaluate a function that wasn't quoted");
  list1 = {"key1": "value1", "key2": "value2"};
  evaluatedList1 = Eval(list1);
  assert(list1["key1"] == evaluatedList1["key1"], "Failed to return the evaluated associative list unchanged");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('Eval(T)', "not implemented/defined for a Topology"), "Failed error checking for trying to take a Evalonential of a topology");
  assert (runCommandWithSoftErrors ('Eval (TT)', "not implemented/defined for a Tree"), "Failed error checking for trying to take a Evalonential of a tree");


  testResult = 1;

  return testResult;
}
