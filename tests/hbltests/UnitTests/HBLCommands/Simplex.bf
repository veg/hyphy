ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Simplex";
}

function runTest ()
{

    testResult = 0;

    m = {
        {                 1,                 2,                 3,                 4,                 0,                 0,                 0}
        {                 0,                 3,                 2,                 1,                 1,                 0,                10}
        {                 0,                 2,                 5,                 3,                 0,                 1,                15}
    };

    result = {{                 1,                 0,                 0,                -0,                 0,                 0}};

    assert(Simplex(m)==result,"Simplex failed");

    testResult = 1;
    return testResult;
}
