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
    assert (Log(export_string) == 3059788636, "Failed checksum test when exporting a likelihood function object");
    ExecuteCommands (export_string);
    Export (export_string, MG94customModel);
    assert (Log(export_string) == 2553492864, "Failed checksum test when exporting a model object");
    Export (export_string, filteredData);
    assert (Log(export_string) == 520896235, "Failed checksum test when exporting a data filter object");
    
    ExecuteCommands (export_string);
    DataSet roundtrip_filter = ReadFromString (export_string);
    assert (roundtrip_filter.sites == 1320 && roundtrip_filter.species == 8, "Failed a dataset filter roundtrip test");
    
    
    assert (runCommandWithSoftErrors ("Export (seamus/romney, some_object)", "is not a valid variable identifier in call to Export"), "Invalid variable identifier in call to Export.");
    assert (runCommandWithSoftErrors ("Export (seamus_romney, some_object)", " is not a supported type in call to Export"), "Invalid object to export.");

	testResult = 1;
		
	return testResult;
}