function getTestName ()
{
	return "Type";
}		

function getTestedFunctions ()
{
	return {{"_MathObject::Type;*::ObjectClass"}};
}	

function runTest ()
{
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
	assert (Type(Exp(1)) == "Number", "Type of a numeric object");

	M = {{1,2}};
	assert (Type(M) == "Matrix", "Type of a matrix object");

	A = {"0":"A", "A":"0"};
	assert (Type(A) == "AssociativeList", "Type of an associative array object");

	Topology Top = (1,2,3);
	assert (Type(Top) == "Topology", "Type of a topology object");

	Tree 	T = (1,2,3);
	assert (Type(T) == "Tree", "Type of a tree object");
	assert (Type (BranchLength (T,5)) == "Unknown", "Type of an unknown object");
	
	assert (Type("Hello world") == "String", "Type of a string object");
	
	testResult  	   = 1;
	return testResult;
}

/* execution stub */

fprintf    (stdout, "[Running COVERAGE TEST '", getTestName(), "']\n");
result  =  runTest();
if (result)
{
	fprintf (stdout, "[TEST PASSED]\n");
}
else
{
	fprintf (stdout, "[TEST FAILED]\n");
}
