ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "#include";
}		

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::ProcessInclude"}};
}	

function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
    #include "nested/nested2/bar.bf";
    assert (fubar=='FOObar', "Checking for the correct inclusion of files");

    #include"nested/baz.bf";
    assert (nonce > 1, "Checking for the correct inclusion of files where there is no space");


    assert (runCommandWithSoftErrors (' #include "nosuchfile";', "Could not read batch file"), "Failed error checking for an invalid file path");
   
    

	testResult = 1;
		
	return testResult;
}
