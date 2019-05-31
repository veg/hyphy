ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Format";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // TODO: Format a number to a specified number of sigfigs and decimal places (the sigfigs argument doesn't seem to have an effect).
  unformattedX = 51.123456;
  formattedX1 = Format(unformattedX, 3, 3);
  formattedX2 = Format(unformattedX, 4, 3);
  formattedX3 = Format(unformattedX, 2, 4);
  assert(formattedX2 == '51.123', "Failed to correctly format a number (Format(51.123456, 4, 3))");

  // SHOWING THAT THE SIGFIGS ARGUMENT DOESN'T SEEM TO HAVE AN EFFECT)
  //fprintf (stdout, 'formattedX1: ', formattedX1, '\n');
  //fprintf (stdout, 'formattedX2: ', formattedX2, '\n');
  //fprintf (stdout, 'formattedX3: ', formattedX3, '\n');//assert(formattedX1 == '51.12', "Failed to correctly format a number (Format(51.123456, 3, 3))");
  //assert(formattedX3 == '51', "Failed to correctly format a number (Format(51.123456, 2, 4))");

  // Test rounding.
  assert(Format(0.116, 2, 2) == '0.12', "Failed to format with rounding");
  assert(Format(0.114, 2, 2) == '0.11', "Failed to format with rounding");
  assert(Format(0.115, 2, 2) == '0.12', "Failed to format with rounding (rounding up a 5)");

  // None is treated like zero.
  assert(Format(none, 2, 2) == '0.00', "Failed to format 'none' as if it were zero");

  // TODO: Inputing a string the string will be treated as a number (according to docs).
  //assert(Format('0.116', 2, 2 == '0.12', "Failed to format a string"));

  // Format a Topology in an easy to read string (topology, internalNodeNames*, branchLengths*) *: if there's a number, show, if not don't.
  Topology T = ((1:.34,b),(3:200,4),5);
  topoTrueTrue = Format(T, 1, 1);
  topoFalseFalse = Format(T, false, false);
  topoTrueFalse = Format(T, 1, false);
  topoFalseTrue = Format(T, false, 1);
  // TODO: (The below tests were passing before ~1/20/2019) 
  //assert(topoTrueTrue == '((1:0.34,b:-1)Node1:-1,(3:200)Node4:-1,5:-1)', "Failed to format topology with internal node names and branch lengths");
  //assert(topoFalseFalse == '((1,b),(3),5)', "Failed to format topology without internal node names and branch lengths");
  //assert(topoTrueFalse == '((1,b)Node1,(3)Node4,5)', "Failed to format topology with internal node names but without branch lengths");
  //assert(topoFalseTrue == '((1:0.34,b:-1):-1,(3:200):-1,5:-1)', "Failed to format topology without internal node names but with branch lengths");
  // Trees are formatted in the same way as topologies
  Tree TT = ((1:.34,b),(3:200,4),5);
  treeTrueTrue = Format(TT, 1, 1);
  treeFalseFalse = Format(TT, false, false);
  treeTrueFalse = Format(TT, 1, false);
  treeFalseTrue = Format(TT, false, 1);
  //assert(treeTrueTrue == '((1:0.34,b:-1)Node1:-1,(3:200)Node4:-1,5:-1)', "Failed to format tree with internal node names and branch lengths");
  //assert(treeFalseFalse == '((1,b),(3),5)', "Failed to format tree without internal node names and branch lengths");
  //assert(treeTrueFalse == '((1,b)Node1,(3)Node4,5)', "Failed to format tree with internal node names but without branch lengths");
  //assert(treeFalseTrue == '((1:0.34,b:-1):-1,(3:200):-1,5:-1)', "Failed to format tree without internal node names but with branch lengths");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "value1", "key2": "value2"};

  // TODO: the below was passing before (~1/20/2019)
  // assert (runCommandWithSoftErrors ('Format(list1)', "not implemented/defined for a AssociativeList"), "Failed error checking for trying to take a Formating of a topology");


  testResult = 1;

  return testResult;
}
