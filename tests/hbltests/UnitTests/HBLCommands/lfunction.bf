ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "lfunction";
}		

lfunction redefine () {
  thisIsNotRedefined = 0;
  return null;
}

lfunction sum (a,b) {
    return a + b;
}

lfunction addTwoByValue (x) {
  x = x + 2;
  return x;
}

lfunction addTwoByRef (y&) {
  y = y + 2;
  return y;
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;

  // Test to confirm that variables declared within an lfunction are at local scope only (function variables {note the lack of the `l` are at global scope})
  thisIsNotRedefined = 1;
  redefine();
  assert(thisIsNotRedefined, "A variable defined in an lfunction leaked into the global scope");

  testSum = sum(3,4);
  assert(testSum == 7, "Failed to successfully define and execute a sum lfunction");

  testNum1 = 5;
  testNum1Plus2 = addTwoByValue(testNum1);
  assert(testNum1Plus2 - testNum1 == 2, "Failed to successfully define and execute an lfunction using parameter by value");

  testNum2 = 6;
  testNum2Plus2 = addTwoByRef('testNum2');
  assert(testNum2Plus2 == testNum2, "Failed to successfully define and execute an lfunction using parameter by refernce");

  testResult = 1;

  return testResult;
}
