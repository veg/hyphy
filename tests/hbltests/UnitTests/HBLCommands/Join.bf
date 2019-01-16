ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Join";
}

function runTest ()
{

    testResult = 0;

    // Test string functionality
    list={"key":"value", "key2":"value2"};
    assert(Join(",", list)=="value,value2","Joining an associative list failed");

    mat = {{1,2,3}{4,5,6}{7,8,9}};
    assert(Join(",", mat)=="1,2,3,4,5,6,7,8,9","Joining a matrix failed");

    testResult = 1;
    return testResult;
}
