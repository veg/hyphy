ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Inverse";
}

function runTest ()
{

    testResult = 0;

    // Test string functionality
    dracula = "dracula";
    reverse = Inverse(dracula);
    assert(Inverse(dracula)=="alucard","Inverse string failed");

    palindrome = "-amanaplanacanalpanama-";
    assert(Inverse(palindrome)==palindrome,"Inverse string failed");


    // Test matrix functionality
    A = {{1,3,3}{1,4,3}{1,3,4}};
    B = Inverse(A);
    I = {{1,0,0}{0,1,0}{0,0,1}};

    assert(A*B==I,"Inverse matrix failed");

    testResult = 1;
    return testResult;
}
