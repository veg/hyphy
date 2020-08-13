ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName ()
{
	return "DeleteObject";
}

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::HandleDeleteObject"}};
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

  ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc2.nex");
  ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc.nex", {}, "a_prefix");
  ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc.nex");
	testResult  	   = TRUE;



	return 1;
}
