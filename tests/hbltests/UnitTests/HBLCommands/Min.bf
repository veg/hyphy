ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName ()
{
  return "Min";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Compare two numbers.
  assert(Min(1,2) == 1, "Failed to compute minimum of two numbers");
  assert(Min(5,4.999) == 4.999, "Failed to compute minimum of two numbers that are similar in magnitude");
  assert(Min(-1,-2) == -2, "Failed to compute minimum of two negative numbers");
  // Find min of array.
  assert(Min({{1,2}{3,4}},2) == 1, "Failed to compute minimum value of an array");
  // Find min of array and it's index. 
  assert(Min({{1,2}{3,4}},1) == {{1, 0}}, "Failed to compute minimum value and index in an array");
  // Min functions on dictionary; should only select among keys with numeric values
  dict = {"0" : 1, "1" : {{2,3}}, "hai" : {"a" : 5, "b" : 7}, "beavis" : 42};
  assert((Min(dict))["key"] == "0", "Failed to compute minimum value in a dictionary");
  
  // For string, Min returns absolute file paths (absolute filepaths aren't straightforward to test so shown here for documentation):
  // with an empty string it returns a temp directory "/tmp/HYPHY-<6Characters_mixed_upper_and_lower_case>"
  tempDir = Min("",0);
  //fprintf (stdout, 'tempDir: ', tempDir, '\n');
  // With a non-empty string it returns the current working directory pluss the string
  // With second argument of zero it gets current working directory of the hbl script
  workDirZero = Min("test",0);
  //fprintf (stdout, 'workDirZero: ', workDirZero, '\n');
  // With the second argument greater than zero it gets current working directory of the hyphy executable
  workDirOne = Min("tests",1);
  //fprintf (stdout, 'workDirOne: ', workDirOne, '\n');

  //---------------------------------------------------------------------------------------------------------
  // TOPOLOGY
  //---------------------------------------------------------------------------------------------------------
  // Can't seem to get the topology function to work... 
  // getting the folowing error: The left hand side expression does not contain an object reference in the following context: 'Min(Topology  T1=<ERROR HERE>(1,2,3);,2)testResult=1'



  Topology T1 = ((1:0.1, 2:0.2)N12 : 0.5, 3 : 1, 4 : 1);
  
  //fprintf (stdout, Min (T1,2), "\n");
 
  // Example from Docs: http://hyphy.org/w/index.php/Min
  //Min(Topology T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);, 2);

  // Simpler example
  //Min (Topology T1 = (1,2);, 2);

  // Other example. The below example returns the following error:
  /*
  'RPN:T|2|min|*' evaluated with errors. Unconsumed values on the stack
  [2]------------------
  0
  [1]------------------
  0
  */
  //Topology T = ((1,2),(3,4),5);
  //min(T,2);


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Tree TT = ((1,2),(3,4),5);
  Tree TTb = ((1,2),(3,4),5,6);

  assert (runCommandWithSoftErrors ('Min (None,0)', "Attempting to operate on an undefined value"), "Failed error checking for trying to take a minimum with None");
  assert (runCommandWithSoftErrors ('Min (TT, TTb)', "Invalid power argument in call to COT finder"), "Failed error checking for trying to take a minimum of a tree");
  assert (runCommandWithSoftErrors ('Min (1)',  "was called with an incorrect number of arguments"), "Too few arguments error check");


  testResult = 1;

	return testResult;
}
