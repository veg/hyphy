ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "*";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Multiply numbers 
  assert(4*4 == 16, "Failed to multiply two integers (4*4)");
  assert(4.1561*4.4562 == 18.52041282, "Failed to multiply two numbers (4.1561*4.4562)");
  assert(Abs(4.156195*4.456201 - 18.5208403152) < 1e-9, "Failed to multiply two decimals to four decimal places");
  assert(4*2*5 == 40, "Failed to multiply 3 numbers");
  
  
  // String buffering
  
  concatedString = "";
  concatedString * "juxta";
  concatedString * "position";
  
  assert(concatedString == "juxtaposition", "Failed to concatenate two strings");
   
  // Multiply each element in a matrix by a number
  assert({{1,2,3}{4,5,6}}*2 == {{2,4,6}{8,10,12}}, "Failed to multiply each element in a matrix by a number");

  // Standard matrix multiplication
  assert({{2,2,4,2}}*{{1}{3}{5}{7}} == {{42}}, "Failed to perform matrix multiplication");

  // String * number clears string, by allocating an empty buffer
  clearedString = "vuala";
  clearedString*2;
  assert(clearedString == "", "Failed to clear a string");

  // Map a string to a vector
  mappedToVector = "tatatgx"*{{"t", "a", "g"}};
  assert(mappedToVector == {{0,1,0,1,0,2,-1}}, "Failed to map a string to a vector");

  // Join two associative lists (update the first list with all the key: value pairs from the second list)
  // duplicate keys receive updated values 
  list1 = {"key1":1, "key2":2};
  list2 = {"key2":3, "key4":4};
  merged_list = list1*list2;
  list2key3 = list2["key2"];
  list1key3 = list1["key2"];
  
  assert(list2["key2"] == list1["key2"] && (list1 / "key4"), "Failed to join two associated lists");
  assert(merged_list == 3, "Failed to return the length of the jointed list when joining two lists");



  
  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  matrix = {{0.3,0.4}};

  // TODO: T*T; this command casues the script to hang.
  //assert (runCommandWithSoftErrors ('T*T', "is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (&&)");
  // TODO: TT*TT; also causes the script to hang.
  //assert (runCommandWithSoftErrors ('TT*TT', "is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (&&)");
  assert (runCommandWithSoftErrors ('4*"String"', "where 'X' is not a number"), "Failed error checking for trying to compare number&&string");

  /* TODO:
  assert(None*4 == 0, "Failed to multiply a number and none");
  assert(4*None == 0, "Failed to multiply a number and none");
  */
  

  testResult = 1;

  return testResult;
}
