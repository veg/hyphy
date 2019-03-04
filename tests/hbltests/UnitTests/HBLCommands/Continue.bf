ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Continue";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Continue is a control flow tool that may only exist during a while or for loop. It continues with the next cycle of the nearest enclosing loop.
  // The below test is designed to check that continue:
  //    works in both for and while loops 
  //    that it exits when encountered
  //    and that it only exits the nearest enclosing loop.

  j = 0;
  x = 0;
  for (y=0; y<=10; y=y+1) {
    for (k = 0;k <= 10; k=k+1) {   
      if (k >= 2)
      {   
          continue;
      }    
      
      j = k;
    }
    x = x+1;
  }
 
  assert(j == 1, "`continue` behaved unexpectedly in a for loop");
  assert(x == 11, "`continue` behaved unexpectedly and seems to have exited more than just the nearest enclosing loop");

  j = 0;
  k = 0;
  while (k<10) {   
    k = k+1;
    if (k >= 2)
    {   

        continue;
    }    
    
    j = k;
  }

  assert(j == 1, "`continue` behaved unexpectedly in a while loop");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('continue;', "'continue;' only makes sense in the context of a loop."), "Failed error checking for trying to execute `continue` outside a loop.");
  assert (runCommandWithSoftErrors ('if(1==1){continue;}', "'continue;' only makes sense in the context of a loop."), "Failed error check for trying to execute `continue` inside an if statment.");
  
  testResult = 1;

  return testResult;
}
