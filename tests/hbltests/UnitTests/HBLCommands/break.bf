
ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "break";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // break is a control flow tool that may only exist during a while or for loop. It terminates nearest enclosing loop.
  
  
  // The below test is designed to check that break:
  //    works in both for and while loops 
  //    that it exits when encountered
  //    and that it only terminates the nearest enclosing loop.
  x = 0;
  y = 10;
  z = 100;

  for (i=0;i<2;i+=1) {
    for (j=0;j<2;j+=1) {
      y+=1;
      k = 0;
      while (k<2) {
        k+=1;
        z+=1;
        break;
      }
      // NOTE: No break statment here.
    }
    break;
    x+=1;
  } 

  assert(x == 0, "`break` behaved unexpectedly");
  assert(y == 12, "`break` behaved unexpectedly");
  assert(z == 102, "`break` behaved unexpectedly");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('break;', "'break;' only makes sense in the context of a loop."), "Failed error checking for trying to execute `break` outside a loop.");
  assert (runCommandWithSoftErrors ('if(1==1){break;}', "'break;' only makes sense in the context of a loop."), "Failed error check for trying to execute `break` inside an if statment.");


  testResult = 1;
  

  return testResult;
}

