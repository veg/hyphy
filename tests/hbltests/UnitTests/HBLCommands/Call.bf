ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Call";
}		

function addTwo(a) {
  return a+2;
}

function sumThreeArgs(a,b,c) {
  return a+b+c;
}

function nonNumericArguments (a,b) {
  return (+a) + (+b);
}

function noArguments () {
  return TRUE;
}

function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  five = Call('addTwo', 3);
  assert(five == 5, "Failed to `Call` a function with one argument as expected");

  six = Call('sumThreeArgs', 1, 2, 3);
  assert(six == 6, "Failed to `Call` a function with multiple arguments as expected");

  forty_two = Call('nonNumericArguments', {{10,9,8,7}}, {"0": 6, "10" : 2});
  assert(forty_two == 42, "Failed to `Call` a function with multiple non-numeric arguments as expected");

  assert(Call ("noArguments") == TRUE, "Failed to `Call` a function without arguments");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  
  assert (runCommandWithSoftErrors ('Call(1);', "Operation 'Call' is not implemented/defined for a Number"), "Failed error checking for executing `Call` on a number");
  assert (runCommandWithSoftErrors ('Call(list1);', "Operation 'Call' is not implemented/defined for a AssociativeList"), "Failed error checking for executing `Call` on an associative list");
  assert (runCommandWithSoftErrors ('Call(T1);', "Operation 'Call' is not implemented/defined for a Topology"), "Failed error checking for executing `Call` on a topology");
  assert (runCommandWithSoftErrors ('Call(TT1);', "Operation 'Call' is not implemented/defined for a Tree"), "Failed error checking for executing `Call` on a tree");
  
  // TODO: When `Call(matrix1);` is run on it's own it has the error listed below... when run in the runCommandWithSoftErrors it has this error, `'Call(matrix1)' evaluated with errors. in call to Call(matrix1);`.
  //assert (runCommandWithSoftErrors ('Call(matrix1);', "Operation 'Call' is not implemented/defined for a Matrix"), "Failed error checking for executing `Call` on a matrix");
  
  testResult = 1;

  return testResult;
}
