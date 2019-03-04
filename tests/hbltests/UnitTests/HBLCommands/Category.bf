ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Category";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;

  category category1 = (4, EQUAL, MEAN, GammaDist(_x_,shapeParameter,shapeParameter), CGammaDist(_x_,shapeParameter,shapeParameter), 0 , 
                        1e25,CGammaDist(_x_,shapeParameter+1,shapeParameter));

  GetInformation(cat1Info, category1);
  fprintf (stdout, 'cat1Info: ', cat1Info, '\n');
  for (rate=0; rate < 4; rate+=10) {
    assert(cat1Info[1][i] == 0.25, "Failed to get equal probability of rates (0.25) for a simple Category variable");
  }

  testResult = 1;

  return testResult;
}
