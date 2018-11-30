ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName ()
{
  return "Log";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Natural logarithm of number.
  assert(Log(1) == 0);
  assert(Log(2.718281828459) > 0.9999999999);
  assert(Log(2.718281828459) < 1.0000000001);
  // Adler-32 checksum of a string.
  assert(Log("hyphy"), 109707827);
  // Natural logarithm of array.
  assert(Log({{1,2}{3,4}}) == {{0,0.6931471805599453}{1.09861228866811,1.386294361119891}});


  testResult = 1;

	return testResult;
}
