ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "FindRoot";
}		

function getTestedFunctions (){
	return {{"_ElementaryCommand::HandleFindRootOrIntegrate","_Formula::Newton","_Formula::Brent"}};
}	

function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
    
    FindRoot (root, x^3-125, x, 0, 10);    
    assert (root == 5, "Failed to solve a differentiable equation");

    FindRoot (root, Abs (x) - 5, x, 0, 10);    
    assert (root == 5, "Failed to solve a non-differentiable equation");

    FindRoot (root, Abs (x) - 5, x, -10, 0);    
    assert (root == -5, "Failed to solve a non-differentiable equation");
    
    y = 3;
    FindRoot (root, Exp (x-y) - 1, x, -10, 10);    
    assert (root == y, "Failed to solve an equation with more than one variable");
    
    
    FindRoot (root, x^3-125, x, y*2, 20);    
    assert (root == y*2, "Did not set root value to left bound when no roots available");

    // error checks here
    

    assert (runCommandWithSoftErrors ('FindRoot (not/id, x, x, -1, 1);', "is not a valid variable identifier"), "Failed error checking for an invalid receptacle name");
    assert (runCommandWithSoftErrors ('FindRoot (root, x#2, x, -1, 1);', "Failed to parse "), "Failed error checking for an invalid expression");
    assert (runCommandWithSoftErrors ('FindRoot (root, x-1/2, y, -1, 1);', "does not depend on the variable "), "Failed error checking expression that does not depend on the target variable");
    assert (runCommandWithSoftErrors ('FindRoot (root, x-1/2, x, 1, -1);', "is not a valid interval"), "Failed error checking invalid interval");
    assert (runCommandWithSoftErrors ('FindRoot (root, x-1/2, x, x#2, 1);', "Unexpected symbol"), "Failed error checking invalid interval specification expression (lb)");
    assert (runCommandWithSoftErrors ('FindRoot (root, x-1/2, x, 1, {{1,2}});', "was expected to be a numerical argument"), "Failed error checking invalid interval specification expression (ub)");


	testResult = 1;
		
	return testResult;
}
