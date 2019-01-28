ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "&&";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Check if two numbers are not zero (TODO: really is checking that both numbers are >=1)
  assert(29&&0 == 0, "Failed to return a 0 for one non-zero and one zero (29&&0)");
  oneNonZero = 0&&5;
  oneLessThanOne = 29&&0.9;
  twoNonZero= 29&&5;
  assert(oneNonZero == 0, "Failed to return a 0 for one non-zero and one zero (0&&5)");
  assert(twoNonZero == 1, "Failed to return a 1 for two non-zero numbers (29&&5)");
  assert(oneLessThanOne == 0, "Failed to return a o for one no-zero and one number less than 1 (29&&0.9)");


  // Modify strings based on the modifier key used like: `stringToModify&&modifierkey`.
  // 1. to UPPER
  upperTest = "hyphy"&&1;
  assert(upperTest == "HYPHY", "Failed to convert string to uppercase");
  // 5. HTML escape
  htmlEscapeTest = "hyphy's"&&5;
  assert(htmlEscapeTest == "hyphy&apos;s", "Failed to convert special characters to HTML format");
  // >6. to lowercase
  lowerTest = "HyPhy"&&7;
  assert(lowerTest == "hyphy", "Failed to convert string to lowercase");

  // TODO: 2 and 3 don't seem to work... see below
  // 2. Escape special characgters with a '\'
  escapeTest = "hyphy\\'"&&2;
  fprintf (stdout, "escpaeTest: ", escpaeTest, "\n");
  escapeTest2 = "hyphy\'"&&2;
  fprintf (stdout, "secapeTest2: ", secapeTest2, "\n");
  // 3. Escape '\n', '\t', '"', ad '\\' with a '\'
  // the `string&&3` command causes the script to hang
  //escapeNewLine = "hyphy\\n"&&3;
  //fprintf (stdout, "escapeNewLine: ", escapeNewLine, "\n");
    

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  matrix = {{0.3,0.4}};

  assert (runCommandWithSoftErrors ('T&&T', "is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (&&)");
  assert (runCommandWithSoftErrors ('TT&&TT', "is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (&&)");
  assert (runCommandWithSoftErrors ('4&&"String"', "where 'X' is not a number"), "Failed error checking for trying to compare number&&string");

  // TODO: Fails error checking for matrices. Not sure why.
  //matrix&&matrix; //Running this commad returns "Error: Invalid baseline p-value (must be in (0,1))"
  //None of the tests below pass
  //assert (runCommandWithSoftErrros ('matrix&&matrix', "Invalid baseline p-value"), "Failed error checking for matrices");
  //assert (runCommandWithSoftErrros ('matrix&&matrix', "evaluated with errors"), "Failed error checking for matrices");
  //assert (runCommandWithSoftErrros ('matrix&&matrix', "Unconsumed values on the stack"), "Failed error checking for matrices");


  testResult = 1;

  return testResult;
}
