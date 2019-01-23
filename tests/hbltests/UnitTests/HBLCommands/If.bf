ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "If";
}		

function runTest () {
	testResult = 0;

  value = 0;
  if(1) {
    value = 1;
  } else {
    value = 2;
  }
  assert(value == 1, "If true should trigger the following block");

  value = 0;
  if(0) {
    value = 1;
  } else {
    value = 2;
  }
  assert(value == 2, "If false should trigger an else block");

  value = 0;
  if("") {
    value = 1;
  }
  assert(value != 1, "If should not trigger with an empty string");

  value = 0;
  if("a") {
    value = 1;
  }
  assert(value == 1, "If should trigger with a non-empty string");


  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  assert(runCommandWithSoftErrors('if({}){}', "did not evaluate to a number, a string, or a null"), "If should not run for an AssociativeList");
  assert(runCommandWithSoftErrors('if(T){}', "did not evaluate to a number, a string, or a null"), "If should not run for a Topology");
  assert(runCommandWithSoftErrors('if(TT){}', "did not evaluate to a number, a string, or a null"), "If should not run for a Tree");
  testResult = 1;

  return testResult;
}
