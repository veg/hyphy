ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "==";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // The "double equals" (`==`) tests for equality, returning 1 if the entities are equal and 0 if they are not.

  // To avoid using the `==` as part of the testing itself I will use >0.5 and <0.5 to test for true or false
  Topology T1 = ((1,2),(3,4),5);
  Topology T2 = ((1,2),(3,4),5);
  Topology T3 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  Tree TT2 = ((1,2),(3,4),5);
  Tree TT3 = ((1,4),(3,4),5);
  matrix1 = {{0.3,0.4}};
  matrix2 = {{0.3,0.4}};
  matrix3 = {{0.4,0.4}};

  assert((1==1)>0.5, "The double equals operator didn't work as expected for 1==1");
  assert((1==0)<0.5, "The double equals operator didn't work as expected for 1==0");
  assert((matrix1==matrix1)>0.5, "The double equals operator didn't work as expected for identical matrices (the actual same object)");
  assert((matrix1==matrix2)>0.5, "The double equals operator didn't work as expected for identical matrices (different objects same values)");
  assert((matrix1==matrix3)<0.5, "The double equals operator didn't work as expected for different matrices");

  //TODO: can't get comparing trees or topologies to work, getting seg faults.
  //assert((T1==T2)>0.5, "The double equals operator didn't work as expected for identical topologies");
  //assert((T1==T3)<0.5, "The double equals operator didn't work as expected for different topologies");
  //assert((TT1==TT2)>0.5, "The double equals operator didn't work as expected for identical trees");
  //assert((TT1==TT3)<0.5, "The double equals operator didn't work as expected for different trees");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  
  // TODO: getting unconsumed values on stack errors from the below:
  //assert (runCommandWithSoftErros ('list1==list1;', "Operation '==' is not implemented/defined for a AssociativeList"), "Failed error checking for comparing associative lists");
  
  testResult = 1;

  return testResult;
}
