ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "GetURL";
}		

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::HandleGetURL", "CheckReceptacleCommandID"}};
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

    GetURL (url_data, "http://www.hyphy.org");
    assert (Abs(url_data) > 0, "Expected to retrieve non-empty data from http://www.hyphy.org");
    
    GetURL (PATH_TO_CURRENT_BF + "tmp" + DIRECTORY_SEPARATOR + "GetURL.txt", "http://www.google.com", SAVE_TO_FILE);
    assert (Abs(url_data) > 0, "Expected to retrieve non-empty data from http://www.hyphy.org");
    
    assert (runCommandWithSoftErrors ("GetURL (url_data/j, \"http://www.google.com\")", " is not a valid variable identifier"), "Failed error checking for an invalid receptacle");
    assert (runCommandWithSoftErrors ("GetURL (url_data, \"http://www.thereisnowaythisurlexists.org\")", "Could not fetch "), "Failed error checking for a URL fetch error");
    assert (runCommandWithSoftErrors ("GetURL (PATH_TO_CURRENT_BF + \"tmp\" + DIRECTORY_SEPARATOR + \"GetURL.txt\", \"http://www.thereisnowaythisurlexists.org\", SAVE_TO_FILE)", "Could not fetch "), "Failed error checking for a URL fetch error, writing to file");
    assert (runCommandWithSoftErrors ("GetURL (url_data, \"http://www.hyphy.org\", FAILED_FLAG)", "Unknown action flag "), "Failed error checking for an invalid action flag");

	testResult = 1;
		
	return testResult;
}
