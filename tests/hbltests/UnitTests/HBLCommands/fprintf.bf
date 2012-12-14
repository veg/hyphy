ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "fprintf";
}		

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::HandleFprintf"}};
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

    fprintf (stdout, "\nTesting basic fprintf functionality. 2+2 = ", 2+2, "\n");
    
    _test_file_name = "tmp/fprintf.txt";
    
    
    fprintf (_test_file_name, CLEAR_FILE, KEEP_OPEN, "A string", {{1,2}}, {"a":"bcd"});
    fprintf (_test_file_name, CLEAR_FILE, KEEP_OPEN, "A string", {{1,2}}, {"a":"bcd"});
    fprintf (_test_file_name, CLOSE_FILE);
    fscanf (_test_file_name, "String,NMatrix,Raw", _s, _m, _a);
    assert (_s == "A string", "fprintf-fscanf roundtrip for failed for a string");
    assert (Abs (m - {{1,2}}) < 1e-12, "fprintf-fscanf roundtrip for failed for a matrix");
    fprintf (_test_file_name, CLEAR_FILE,PRINT_SELF, LIST_ALL_VARIABLES);

    _m = {{1,2}};
    _s = "_m";
    fprintf (TEMP_FILE_NAME, "This will not live for very long", _s__, "\n");
    
    GLOBAL_FPRINTF_REDIRECT = "/dev/null";
    fprintf (stdout, "Should not be displayed");

    GLOBAL_FPRINTF_REDIRECT = "tmp/redirect.txt";
    fprintf (stdout, CLEAR_FILE, "Gone to tmp/redirect.txt");
    GLOBAL_FPRINTF_REDIRECT = "";
    
    assert (runCommandWithSoftErrors ("fprintf (\"tmp/tmp/tmp/can't open this\", \"No go\");", "Could not create/open output file at path"), "Failed error checking for an invalid path specification");
    assert (runCommandWithSoftErrors ("fprintf (stdout, a=b);", "Argument 1 is not a simple expression"), "Failed error checking for trying to print an assignment");
    assert (runCommandWithSoftErrors ("fprintf (stdout, /a);", "Bad binary operator placement "), "Failed error checking for trying to print a borked expression");
    
    
    /*
    GetURL (url_data, "http://www.hyphy.org");
    assert (Abs(url_data) > 0, "Expected to retrieve non-empty data from http://www.hyphy.org");
    
    GetURL (PATH_TO_CURRENT_BF + "tmp" + DIRECTORY_SEPARATOR + "GetURL.txt", "http://www.google.com", SAVE_TO_FILE);
    assert (Abs(url_data) > 0, "Expected to retrieve non-empty data from http://www.hyphy.org");
    
    assert (runCommandWithSoftErrors ("GetURL (url_data/j, \"http://www.google.com\")", " is not a valid variable identifier"), "Failed error checking for an invalid receptacle");
    assert (runCommandWithSoftErrors ("GetURL (url_data, \"http://www.thereisnowaythisurlexists.org\")", "Could not fetch "), "Failed error checking for a URL fetch error");
    assert (runCommandWithSoftErrors ("GetURL (PATH_TO_CURRENT_BF + \"tmp\" + DIRECTORY_SEPARATOR + \"GetURL.txt\", \"http://www.thereisnowaythisurlexists.org\", SAVE_TO_FILE)", "Could not fetch "), "Failed error checking for a URL fetch error, writing to file");
    assert (runCommandWithSoftErrors ("GetURL (url_data, \"http://www.hyphy.org\", FAILED_FLAG)", "Unknown action flag "), "Failed error checking for an invalid action flag");
    */
    
	testResult = 1;
		
	return testResult;
}