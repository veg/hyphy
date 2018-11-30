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
  assert(Min(1,2) == 1);
  assert(Min(5,4.999) == 4.999);
  assert(Min(-1,-2) == -2);
  // Find min of array.
  assert(Min({{1,2}{3,4}},2) == 1);
  // Find min of array and it's index. 
  assert(Min({{1,2}{3,4}},1) == {{1, 0}});


  //---------------------------------------------------------------------------------------------------------
  // TOPOLOGY
  //---------------------------------------------------------------------------------------------------------
  // Can't seem to get the topology function to work... 
  // getting the folowing error: The left hand side expression does not contain an object reference in the following context: 'Min(Topology  T1=<ERROR HERE>(1,2,3);,2)testResult=1'
  // Example from Docs: http://hyphy.org/w/index.php/Min
  //Min(Topology T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);, 2);
  // Simpler example
  //Min (Topology T1 = (1,2);, 2);


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  /*
  Trying to figure out how to use "runCommandWithSoftErrors...
  Using the command: Min({{1,2}{3,4}}) as an example.
  This command is missing the second argument and thus returns the following error to stdout when called:
  "Operation 'Min' was called with an incorrect number of arguments (0) for Matrix"
  */

  //Command that produces the error
  //Min({{1,2}{3,4}});
  
  /*
   Execute the runCommandWithSoftErrors on it's own:
     prints:
        Expected an error matching 'Operation 'Min' was called with an incorrect number of arguments (0) for Matrix', while executing 'Min({{1,2}{3,4}});'.
        Had error
        'Min({
        {1, 2}
        {3, 4}
        })' evaluated with errors. in call to Min({{1,2}{3,4}});
  */
  //runCommandWithSoftErrors("Min({{1,2}{3,4}})", "Operation 'Min' was called with an incorrect number of arguments (0) for Matrix");

  // Try to wrap it in an assert.
  //assert( runCommandWithSoftErrors("Min({{1,2}{3,4}})", "Operation 'Min' was called with an incorrect number of arguments (0) for Matrix"), "What Text goes here?");


  testResult = 1;

	return testResult;
}
