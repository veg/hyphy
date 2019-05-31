ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "&";
}		

lfunction f1 () {
    x = 0;
    return &x;
}

lfunction f2 () {
    x = 1;
    return &x;

}

function runTest () {
  ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
  testResult = FALSE;
    // The ampersand in front of a variable means to take the address of a variable including the namespace.

  namespace test {
      bar = "baz";
      assert(&bar == "test.bar", "Failed to call a variable address (using `&`) as expected from within a namespace");
      assert(bar == "baz", "Failed to call a variable as expected within a namespace");
  }

  assert(&bar == "bar", "Failed to call a variable address (using `&`) as expected from outside the namespace where the variable was defined");
  assert(bar == 0, "Failed to call a call a variable as expected from outside the namespace where the variable was devined");

  foo = "phyHi";
  assert(&foo == "foo", "Failed to call a variable address (using `&`) as expected from the global namespace");
  assert(foo == "phyHi", "Failed to call a call a variable as expected from the global namespace");

  assert(f1 () !=  f2(), "Clash in local name spaces detected");

  testResult = TRUE;

  return testResult;
}
