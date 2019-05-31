ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Return";
}		


function testSecondReturn (x) {
    x = x + 2;
    return "second return shouldn't execute";
    return "the second return executed?!?!";
}

function testIfInFunction (x) {
  if (x == 'test return in if') {
    return 'worked';
  }
  x = x + 2;
  return x;
}

function testLoop (x) {
  while (x < 6) {
    x + 1;
    return 'x less than 6';
  }
  return 'x greater than 6';
}

function testEmptyReturn (x) {
  return ;
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;

  assert(testSecondReturn(4) == "second return shouldn't execute", "Failed to only execute the first return statment in a function with multiple returns");

  assert(testIfInFunction(4) == 6, "Failed to execute the return after an if statment");
  assert(testIfInFunction('test return in if') == 'worked', "Failed to execute the return within an if statment inside a function");

  assert(testLoop(2) == 'x less than 6', "Failed to execute the return in a loop");
  assert(testLoop(7) == 'x greater than 6', "Failed to execute the return after a loop when the loop was never entered");

  assert(testEmptyReturn(2) == null, "failed to return null as an empty return");

  testResult = 1;

  return testResult;
}
