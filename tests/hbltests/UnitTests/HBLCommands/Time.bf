ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Time";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Time(x). If x is a real number, If value >= 1, then return unix timestamp. If value==0, return process' CPU time.
  // TODO: the documentation should be updated to reflect that value doesn't need to equal 0 for the second (cpu time) functionality, it just needs to be less than 1.
  assert(Time(1) > 1546300800, "Failed to provide a current unix timestamp after 01/01/2019");
  assert(Time(1) < 1893456000, "Failed to provide a current unix timestamp before 01/01/2040 (I'll gladly update this test in 21 years)");

  cpuTime1 = Time(0);
  for (x=0;x<100000;x+=1) {
      y = 100;
      z = y^4;
  }
  cpuTime2 = Time(0);
  assert(cpuTime2 - cpuTime1 > 0.01, "Failed to time a process correctly (100^4 calculated 100,000 times took less than a hundredth of a second)");
  assert(cpuTime2 - cpuTime1 < 2, "Failed to time a process correctly (100^4 calculated 100,000 times took more than two seconds");

  // Calling time with undefined or none results in the same behaviour as Time(0).
  cpuTime3 = Time(none);
  cpuTime4 = Time(undefinedVariable);
  assert(cpuTime3 < 10, "Failed to execute Time(0) when running Time(none)");
  assert(cpuTime3 > cpuTime1, "Failed to execute Time(0) when running Time(none)");
  assert(cpuTime4 > cpuTime1, "Failed to execute Time(0) when running Time(undefinedVariable)");
   

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------

  exampleString = "abcd";
  Tree exampleTree = ((1,2),(3,4),5);
  Topology exampleTopology = ((1,2),(3,4),5);
  exampleList = {"key1": "val1", "key2": "val2"};

  assert (runCommandWithSoftErrors ('Time(exampleString)', "Operation 'Time' is not implemented/defined for a String"), "Failed error checking for trying to run Time(string)");
  assert (runCommandWithSoftErrors ('Time(exampleList)', "Operation 'Time' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to run Time(AssociativeList)");
  assert (runCommandWithSoftErrors ('Time(exampleTree)', "Operation 'Time' is not implemented/defined for a Tree"), "Failed error checking for trying to run Time(Tree)");
  assert (runCommandWithSoftErrors ('Time(exampleTopology)', "Operation 'Time' is not implemented/defined for a Topology"), "Failed error checking for trying to run Time(Topology)");


  testResult = 1;
  

  return testResult;
}

