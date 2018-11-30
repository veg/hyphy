ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName ()
{
  return "Max";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Compare two numbers.
  assert(Max(1,2) == 2);
  assert(Max(4,4.001) == 4.001);
  assert(Max(-1,-2) == -1);
  // Find max of array.
  assert(Max({{1,2}{3,4}},2) == 4);
  // Find max of array and it's index. 
  assert(Max({{1,2}{3,4}},1) == {{4, 3}});


  testResult = 1;

	return testResult;
}
