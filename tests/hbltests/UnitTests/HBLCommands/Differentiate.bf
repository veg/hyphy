ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "Differentiate";
}		

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::HandleDifferentiate"}};
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
    
    assert (runCommandWithSoftErrors ("Differentiate (2 invalid, x^2,x)", "is not a valid variable identifier in call to Differentiate"), "Invalid variable identifier in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, x^2,x,\"beavis\")", "The number of times to differentiate must be a non-negative integer in call to Differentiate"), "Invalid order option in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, x^2,x,z#4)", "<ERROR HERE>"), "Unparseable order option in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, Time(x),x)","Differentiation of .+ failed"), "Unparseable order option in call to Differentiate.");

    Differentiate (dfx, x^2,x);
    x = -2; 
    assert (dfx == -4, "Checking (x^2)' == 2x derivative");
    Differentiate (dfx, x^3,x,2);
    assert (Abs(dfx+12)<1e-10, "Checking (x^3)'' == 6x derivative");

	testResult = 1;
		
	return testResult;
}