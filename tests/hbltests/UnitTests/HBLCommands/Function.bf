ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Function";
}		

lfunction redefine () {
  thisIsNotRedefined = 1;
  return null;
}

function sum (a,b) {
    return a + b;
}

function addTwoByValue (x) {
  x = x + 2;
  return x;
}

function addTwoByRef (y&) {
  y = y + 2;
  return y;
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;

  /* TODO: It looks like `function` is only setting variables at local scope, not global scope...
  // Test to confirm that variables declared within a function are at global scope (lfunction variables are at local scope)
  thisIsNotRedefined = 0;
  redefine();
  assert(thisIsNotRedefined, "A variable defined in an function is not available in the global scope");
  */

  testSum = sum(3,4);
  assert(testSum == 7, "Failed to successfully define and execute a sum function");

  testNum1 = 5;
  testNum1Plus2 = addTwoByValue(testNum1);
  assert(testNum1Plus2 - testNum1 == 2, "Failed to successfully define and execute a function using parameter by value");

  testNum2 = 6;
  testNum2Plus2 = addTwoByRef('testNum2');
  assert(testNum2Plus2 == testNum2, "Failed to successfully define and execute a function using parameter by refernce");

  //TODO: Overloading a function causes an 'Unconsumed values on the stack' error:
  //testOverload = sum(1,2,3);

  testResult = 1;

  return testResult;
}
