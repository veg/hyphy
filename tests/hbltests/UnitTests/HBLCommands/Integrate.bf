ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "Integrate";
}		

function getTestedFunctions (){
	return {{"_ElementaryCommand::HandleFindRootOrIntegrate","_Formula::Newton","_Formula::Brent"}};
}	

function runTest () {

    
     

	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
	Beta (0,0);
    
    Integrate (int, x, x, 0, 1);    
    assert (Abs(int - 0.5) < 1e-8, "Failed to integrate x over [0,1]");

    Integrate (int, 1/x^3,x, 1, 1e10);  /* this is, of course -x^(-2) / 2 evaluated between infty and 1, i.e. 0.5 */
    
    // TODO: the below test failes.
    // assert (Abs (int - 0.5) < 1e-4, "Failed integrating a converging function over [1, infty]: " + int);

    for (k = 0; k < 10; k += 1) {
        x = Random (1, 10);
        y = Random (1, 10);
        
        Integrate (int, t^(x-1)*(1-t)^(y-1),t, 1e-10, 1);   // should yield the Beta function B(x,y)
    
        // TODO: the below test fails.
        //assert (Abs ((int - Beta (x,y))/int) < 1e-4, "Failed to integrate the Beta function :" + int + " " + Beta (x,y) + " " + x + " " + y);
    }
    
    for (k = 0; k < 10; k += 1) {
        x = Random (1, 10);
        
        Integrate (int, t^(x-1)*Exp (-t),t, 0, 1e10);   // should yield the Beta function B(x,y)
    
        // TODO: the below test fails.
        //assert (Abs ((int - Gamma (x))/int) < 1e-4, "Failed to integrate the Gamma function :" + int + " " + Gamma (x) + " " + x);
    }


    Integrate (int, 1 ,t, 0, 2);  /* check constant integration */
    assert (Abs (int-2) < 1e-4, "Failed integrating a constant " + int);

    x = 2;
    Integrate (int, x ,t, 0, 2);  /* check constant integration */
    assert (Abs (int-4) < 1e-4, "Failed integrating a function that does not depend on the variable being integrated over " + int);

    // error checks here
    

    assert (runCommandWithSoftErrors ('Integrate (not/id, x, x, 0, 1);', "is not a valid variable identifier"), "Failed error checking for an invalid receptacle name");
    assert (runCommandWithSoftErrors ('Integrate (int, x#2, x, -1, 1);', "Failed to parse "), "Failed error checking for an invalid expression");
    assert (runCommandWithSoftErrors ('Integrate (int, x-1/2, x, 1, -1);', "is not a valid interval"), "Failed error checking invalid interval");
    assert (runCommandWithSoftErrors ('Integrate (int, x-1/2, x, x#2, 1);', "Unexpected symbol"), "Failed error checking invalid interval specification expression (lb)");
    assert (runCommandWithSoftErrors ('Integrate (int, x-1/2, x, 1, {{1,2}});', "was expected to be a numerical argument"), "Failed error checking invalid interval specification expression (ub)");
  

	testResult = 1;
		
	return testResult;
}
