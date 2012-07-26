ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "Export";
}		

function getTestedFunctions ()
{
	return {{"_ElementaryCommand::HandleExport"}};
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
    
    ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc.nex");
    
    Export (export_string, lf);
    DeleteObject (lf);
    ExecuteCommands (export_string);
    assert (trapAllErrors("Export (export_string, MG94customModel);") == 1, "Roundtrip likelihood function export (model)");
    assert (trapAllErrors("Export (export_string, filteredData);") == 1, "Roundtrip likelihood function export (data)");
    
    DataSet roundtrip_filter = ReadFromString (export_string);
    assert (roundtrip_filter.sites == 1320 && roundtrip_filter.species == 8, "Failed a dataset filter roundtrip test");
    
    
    assert (runCommandWithSoftErrors ("Export (seamus/romney, some_object)", "is not a valid variable identifier in call to Export"), "Invalid variable identifier in call to Export.");
    assert (runCommandWithSoftErrors ("Export (seamus_romney, some_object)", " is not a supported type in call to Export"), "Invalid object to export.");

	testResult = 1;
		
	return testResult;
}