ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "RequireVersion";
}		

function runTest () {
  ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
  testResult = TRUE;
  RequireVersion ("1.00");
  RequireVersion ("2.4.0");

  assert(runCommandWithSoftErrors ("RequireVersion ('9.99')", "Current script requires at least version 9.99 of HyPhy"), "Failed to return the expected error message when requiring a hyphy version greater than the current");

  testResult = 1;

  return testResult;
}
