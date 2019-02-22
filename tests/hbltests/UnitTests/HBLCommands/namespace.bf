x = "tu-ti-tu";

namespace foo {
    namespace bar {
        x = "foo";
        namespace baz {
            x = "foo";
            
            function sum (a,b) {
                return a + b;
            }
            
            function inner_wrapper () {
                return echo();
            }
        }
    }
    
    function echo () {
        return "top level";
    }
    
    function wrapper (a,b) {
        return bar.baz.sum (a,b);
    }
    
    x = "foo";
    
}
ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName () {
	return "namespace";
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
    testResult = 0;
    
	assert (foo.x == "foo", "Failed scoping outer namespace (had '" + foo.x + "'" + ")");
	assert (foo.bar.x == "foo", "Failed scoping second-level namespace (had '" + foo.bar.x + "'" + ")");
	assert (foo.bar.baz.x == "foo", "Failed scoping inner namespace (had '" + foo.bar.baz.x + "'" + ")");
	assert (x != "foo", "Incorrectly initialized global scope");
	assert (foo.wrapper (2,3) == 5, "Incorrect result from a wrapper function");
	assert (foo.bar.baz.sum (2,3) == 5, "Incorrect result from a fully resolved function");
	assert (foo.bar.baz.inner_wrapper () == foo.echo (), "Incorrect result from a nested function calling a higher level function");
	
	namespace foo.bar.baz {
	    assert (sum (2,3) == 5, "Local namespace resolution fail");
	}
	
    assert (runCommandWithSoftErrors ("namespace bad-id {}", "Not a valid function/namespace identifier"), "Failed error checking for an invalid namespace ID");

    testResult = 1;
	return testResult;
}
